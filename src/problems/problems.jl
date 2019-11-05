using Plots

export process!, initialize, apply_boundary_conditions!,
    density,
    velocity,
    pressure,
    temperature,
    decay,
    force,
    initialize,
    InitialValueProblem, viscosity, delta_t,
    lattice_velocity,
    lattice_density,
    lattice_pressure,
    lattice_force,
    lattice_viscosity

abstract type LBMProblem end
abstract type SteadyStateProblem <: LBMProblem end
abstract type TimeDependantProblem <: LBMProblem end

abstract type InitialValueProblem end
abstract type DoubleDistributionProblem <: InitialValueProblem end

has_external_force(problem::InitialValueProblem) = false

function initial_equilibrium(quadrature::Quadrature, problem::InitialValueProblem, x::Float64, y::Float64)
    return equilibrium(
        quadrature,
        lattice_density(quadrature, problem, x, y),
        lattice_velocity(quadrature, problem, x, y),
        lattice_temperature(quadrature, problem, x, y)
    )
end

function initial_condition(q::Quadrature, problem::InitialValueProblem, x::Float64, y::Float64)
    initial_equilibrium(q, problem, x, y)
end

import Base: range
function range(problem::InitialValueProblem)
    Δx = problem.domain_size[1] / problem.NX
    Δy = problem.domain_size[2] / problem.NY

    x_range = range(Δx / 2, problem.domain_size[1] - Δx / 2, length = problem.NX)
    y_range = range(Δy / 2, problem.domain_size[1] - Δy / 2, length = problem.NY)


    return x_range, y_range
end

function initialize(quadrature::Quadrature, problem::InitialValueProblem)
    # force_field = Array{Any}(undef, problem.NX, problem.NY)
    f = Array{Float64}(undef, problem.NX, problem.NY, length(quadrature.weights))

    # NOTE: we have periodic boundaries
    x_range, y_range = range(problem)
    # internal_nodes_index(problem)
    for x_idx in 1:problem.NX, y_idx in 1:problem.NY
        f[x_idx, y_idx, :] = initial_equilibrium(
            quadrature,
            problem,
            x_range[x_idx],
            y_range[y_idx]
        )

        # force_field[x_idx, y_idx, :] = force(problem, x, y)
    end
    force_field = (x_idx, y_idx, t) -> lattice_force(problem, x_idx, y_idx, t)

    τ = quadrature.speed_of_sound_squared * lattice_viscosity(problem) + 0.5
    @show τ

    if has_external_force(problem)
        collision_operator = SRT_Force(τ, force_field)
    else
        collision_operator = SRT(τ)
    end

    # if (
    #     typeof(problem) == LinearizedThermalDiffusion ||
    #     typeof(problem) == LinearizedTransverseShearWave
    # )
    # end

    return f, collision_operator
end

function apply_boundary_conditions!(q::Quadrature, problem::InitialValueProblem, f_in, f_out; time = 0.0)
    nothing
end
function apply_boundary_conditions_before!(q::Quadrature, problem::InitialValueProblem; time = 0.0, f_new, f_old)
    nothing
end
function apply_boundary_conditions_after!(q::Quadrature, problem::InitialValueProblem; time = 0.0, f_new, f_old)
    nothing
end


function force(problem::InitialValueProblem, x_idx::Int64, y_idx::Int64, time::Float64 = 0.0)
    x_range, y_range = range(problem)

    x = x_range[x_idx]
    y = y_range[y_idx]

    return force(problem, x, y, time)
end

function is_steady_state(problem::InitialValueProblem)
    return problem.static
end
function is_time_dependant(problem::InitialValueProblem)
    return ! problem.static
end

# Dimensionless
function viscosity(problem) #::InitialValueProblem)
    return problem.ν * delta_x(problem)^2 / delta_t(problem)
end

function heat_diffusion(problem) #::InitialValueProblem)
    return problem.κ * delta_x(problem)^2 / delta_t(problem)
end

function delta_t(problem::InitialValueProblem)
    return delta_x(problem) * problem.u_max
end

function delta_x(problem::InitialValueProblem)
    return problem.domain_size[1] * (1 / problem.NX)
end
function reynolds(problem::InitialValueProblem)
    return problem.NY * problem.u_max / problem.ν
end

lattice_viscosity(problem) = problem.ν #::InitialValueProblem)
lattice_density(q, problem::InitialValueProblem, x, y, t = 0.0) = density(q, problem, x, y, t)
lattice_velocity(q, problem::InitialValueProblem, x, y, t = 0.0) = problem.u_max * velocity(problem, x, y, t)
lattice_pressure(q, problem::InitialValueProblem, x, y, t = 0.0) = problem.u_max^2 * pressure(q, problem, x, y, t)
lattice_force(problem::InitialValueProblem, x, y, t = 0.0) = problem.u_max * delta_t(problem) * force(problem, x, y, t)
function lattice_temperature(q, problem::InitialValueProblem, x, y, t = 0.0)
    return pressure(q, problem, x, y) / density(q, problem, x, y)


    θ = pressure(q, problem, x, y, t) / density(q, problem, x, y, t)

    return θ * problem.u_max^2 #/ q.speed_of_sound_squared
    return θ * problem.u_max^2 / q.speed_of_sound_squared
        # 1.0 / quadrature.speed_of_sound_squared
end

dimensionless_viscosity(problem) = problem.ν * delta_x(problem)^2 / delta_t(problem)
dimensionless_density(problem::InitialValueProblem, ρ) = ρ
dimensionless_velocity(problem::InitialValueProblem, u) = u / problem.u_max
dimensionless_pressure(q, problem::InitialValueProblem, p) = p# (p - 1.0) / (q.speed_of_sound_squared * problem.u_max^2 )
dimensionless_temperature(q, problem::InitialValueProblem, T) = T #* q.speed_of_sound_squared
dimensionless_force(problem::InitialValueProblem, F) = F / (problem.u_max * delta_t(problem))

include("taylor-green-vortex-decay.jl")
include("decaying-shear-flow.jl")
include("poiseuille.jl")
include("couette-flow.jl")
include("lid-driven-cavity.jl")
include("linear-hydrodynamics-modes.jl")
include("convergence.jl")

# error(::Val{:density}, node, solution) = density(node) - density(solution)

struct LatticeProblem2D
    NX::Int
    NY::Int
    τ::Float64
    u_max::Float64
end

