abstract type FluidFlowProblem end
abstract type SteadyStateProblem <: FluidFlowProblem end
abstract type TimeDependantProblem <: FluidFlowProblem end

abstract type InitialValueProblem <: FluidFlowProblem end
abstract type DoubleDistributionProblem <: FluidFlowProblem end

has_external_force(::FluidFlowProblem) = false
function velocity_gradient(
    ::FluidFlowProblem,
    x::Float64,
    y::Float64,
    timestep::Float64 = 0.0,
)
    return zeros(2, 2)
end

function range(problem::FluidFlowProblem)
    Δx = problem.domain_size[1] / problem.NX
    Δy = problem.domain_size[2] / problem.NY

    x_range = range(Δx / 2, stop = problem.domain_size[1] - Δx / 2, length = problem.NX)
    y_range = range(Δy / 2, stop = problem.domain_size[1] - Δy / 2, length = problem.NY)

    return x_range, y_range
end

boundary_conditions(problem::FluidFlowProblem) = BoundaryCondition[]

function deviatoric_tensor(
    q::Quadrature,
    problem::FluidFlowProblem,
    x::Float64,
    y::Float64,
    time::Float64 = 0.0,
)
    a = velocity_gradient(problem, x, y, time)
    ν = viscosity(problem)

    σ = - ν * [
        2 * a[1, 1] a[1, 2] + a[2, 1]
        a[1, 2] + a[2, 1] 2 * a[2, 2]
    ]

    return σ
end

function pressure_tensor(
    q::Quadrature,
    problem::FluidFlowProblem,
    x::Float64,
    y::Float64,
    time::Float64 = 0.0,
)
    A = problem.A
    B = problem.B

    p = pressure(q, problem, x, y, time) * I
    return p - deviatoric_tensor(q, problem, x, y, time)
end

function force(problem::FluidFlowProblem, x_idx::Int64, y_idx::Int64, time::Float64 = 0.0)
    x_range, y_range = range(problem)

    x = x_range[x_idx]
    y = y_range[y_idx]

    return force(problem, x, y, time)
end

is_steady_state(problem::FluidFlowProblem) = problem.static
is_time_dependant(problem::FluidFlowProblem) = !problem.static

# Dimensionless
viscosity(problem::FluidFlowProblem) = problem.ν * delta_x(problem)^2 / delta_t(problem)
heat_diffusion(problem::FluidFlowProblem) = problem.κ * delta_x(problem)^2 / delta_t(problem)
reynolds(problem::FluidFlowProblem) = problem.NY * problem.u_max / problem.ν

function delta_t(problem::FluidFlowProblem)
    return delta_x(problem) * problem.u_max
end

function delta_x(problem::FluidFlowProblem)
    if (problem.NX > problem.NY)
        return problem.domain_size[1] * (1 / problem.NX)
    end

    return problem.domain_size[2] * (1 / problem.NY)
end

lattice_viscosity(problem) = problem.ν #::FluidFlowProblem)
lattice_density(q, problem::FluidFlowProblem, x, y, t = 0.0) = density(q, problem, x, y, t)
lattice_velocity(q, problem::FluidFlowProblem, x, y, t = 0.0) =
    problem.u_max * velocity(problem, x, y, t)
lattice_pressure(q, problem::FluidFlowProblem, x, y, t = 0.0) =
    problem.u_max^2 * pressure(q, problem, x, y, t)
lattice_force(problem::FluidFlowProblem, x, y, t = 0.0) =
    problem.u_max * delta_t(problem) * force(problem, x, y, t)
lattice_temperature(q, problem::FluidFlowProblem, x, y, t = 0.0) =
    pressure(q, problem, x, y) / density(q, problem, x, y)

dimensionless_viscosity(problem) = problem.ν * delta_x(problem)^2 / delta_t(problem)
dimensionless_density(problem::FluidFlowProblem, ρ) = ρ
dimensionless_velocity(problem::FluidFlowProblem, u) = u / problem.u_max
dimensionless_pressure(q, problem::FluidFlowProblem, p) = p# (p - 1.0) / (q.speed_of_sound_squared * problem.u_max^2 )
dimensionless_temperature(q, problem::FluidFlowProblem, T) = T #* q.speed_of_sound_squared
dimensionless_force(problem::FluidFlowProblem, F) = F / (problem.u_max * delta_t(problem))
dimensionless_stress(problem::FluidFlowProblem, σ) = begin
    # Rescale to dimensionless number (TODO check why problem.u_max)
    factor = 1 / problem.u_max^2
    # factor *= 0.964234137002314

    return σ * factor
end

include("taylor-green-vortex.jl")
include("decaying-shear-flow.jl")
include("poiseuille.jl")
include("couette-flow.jl")
include("lid-driven-cavity.jl")
include("linear-hydrodynamics-modes.jl")
include("second-order-convergence.jl")

# error(::Val{:density}, node, solution) = density(node) - density(solution)
