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

function initial_condition(q::Quadrature, problem::FluidFlowProblem, x::Float64, y::Float64)
    f = equilibrium(
        q,
        lattice_density(q, problem, x, y),
        lattice_velocity(q, problem, x, y),
        lattice_temperature(q, problem, x, y),
    )

    q = q
    ρ = sum(f)
    cs = 1 / q.speed_of_sound_squared
    d_u = -0.5 * cs * problem.u_max * lbm.velocity_gradient(problem, x, y, 0.0)
    τ = q.speed_of_sound_squared * lbm.lattice_viscosity(problem) + 0.5
    N = 2

    u = lattice_velocity(q, problem, x, y)
    T = lattice_temperature(q, problem, x, y)
    a_2_eq = equilibrium_coefficient(Val{2}, q, ρ, u, T)

    τ = problem.ν * cs
    ρT = 1.0
    # TODO Use deviatoric_stress_tensor instead?
    a_bar_2 =
        -(1 + 1 / (2 * τ)) *
        τ *
        ρT *
        problem.u_max^2 *
        lbm.velocity_gradient(problem, x, y, 0.0) - (1 / (2 * τ)) * a_2_eq

    return f
    for f_idx = 1:length(f)
        # f[f_idx] += - (q.weights[f_idx] * ρ * τ / cs) * dot(lbm.hermite(Val{2}, q.abscissae[:, f_idx], q), d_u)
        cs = q.speed_of_sound_squared
        f[f_idx] +=
            q.weights[f_idx] *
            (cs^2 / factorial(2)) *
            dot(a_bar_2, hermite(Val{2}, q.abscissae[:, f_idx], q))
    end
    return f
    initial_equilibrium(q, problem, x, y)
end

function range(problem::FluidFlowProblem)
    Δx = problem.domain_size[1] / problem.NX
    Δy = problem.domain_size[2] / problem.NY

    x_range = range(Δx / 2, stop = problem.domain_size[1] - Δx / 2, length = problem.NX)
    y_range = range(Δy / 2, stop = problem.domain_size[1] - Δy / 2, length = problem.NY)

    return x_range, y_range
end

function initialize(quadrature::Quadrature, problem::FluidFlowProblem, cm = SRT)
    f = Array{Float64}(undef, problem.NX, problem.NY, length(quadrature.weights))

    x_range, y_range = range(problem)
    for x_idx = 1:problem.NX, y_idx = 1:problem.NY
        f[x_idx, y_idx, :] =
            initial_condition(quadrature, problem, x_range[x_idx], y_range[y_idx])
    end

    return f
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

    return ν * [
        2 * a[1, 1] a[1, 2] + a[2, 1]
        a[1, 2] + a[2, 1] 2 * a[2, 2]
    ]
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
    return p + deviatoric_tensor(q, problem, x, y, time)
end

function force(problem::FluidFlowProblem, x_idx::Int64, y_idx::Int64, time::Float64 = 0.0)
    x_range, y_range = range(problem)

    x = x_range[x_idx]
    y = y_range[y_idx]

    return force(problem, x, y, time)
end

function is_steady_state(problem::FluidFlowProblem)
    return problem.static
end
function is_time_dependant(problem::FluidFlowProblem)
    return !problem.static
end

# Dimensionless
function viscosity(problem) #::FluidFlowProblem)
    return problem.ν * delta_x(problem)^2 / delta_t(problem)
end

function heat_diffusion(problem) #::FluidFlowProblem)
    return problem.κ * delta_x(problem)^2 / delta_t(problem)
end

function delta_t(problem::FluidFlowProblem)
    return delta_x(problem) * problem.u_max
end

function delta_x(problem::FluidFlowProblem)
    if (problem.NX > problem.NY)
        return problem.domain_size[1] * (1 / problem.NX)
    end

    return problem.domain_size[2] * (1 / problem.NY)
end
function reynolds(problem::FluidFlowProblem)
    return problem.NY * problem.u_max / problem.ν
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

# error(::Val{:density}, node, solution) = density(node) - density(solution)
