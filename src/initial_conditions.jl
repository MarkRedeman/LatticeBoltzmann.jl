abstract type InitializationStrategy end

# TODO: base this on the existence of a velocity gradient, pressure and velocity
# InitializationStrategy(problem::FluidFlowProblem) = AnalyticalEquilibriumAndOffEquilibrium()
InitializationStrategy(problem::FluidFlowProblem) = AnalyticalEquilibrium()

function initialize(strategy::InitializationStrategy, q::Quadrature, problem::FluidFlowProblem, cm = SRT)
    f = Array{eltype(q.weights)}(undef, problem.NX, problem.NY, length(q.weights))

    x_range, y_range = range(problem)
    for x_idx = 1:problem.NX, y_idx = 1:problem.NY
        f[x_idx, y_idx, :] = initial_condition(
            strategy,
            q,
            problem,
            x_range[x_idx],
            y_range[y_idx]
        )
    end

    return f
end
initial_condition(
    strategy::InitializationStrategy,
    q::Quadrature,
    problem::FluidFlowProblem,
    x::T,
    y::T
) where { T <: Real } = initial_condition(q, problem, x, y)

include("initial_conditions/constant_density.jl")
include("initial_conditions/analytical_equilibrium.jl")
include("initial_conditions/analytical_offequilibrium.jl")
include("initial_conditions/analytical_velocity.jl")
include("initial_conditions/analytical_velocity_stress.jl")

"""
f₀(x) = f_eq(ρ₀, u₀, p₀)
ρ₀ = 1.0
u₀ = 0.0
p₀ = 1.0
"""
struct ZeroVelocityInitialCondition <: InitializationStrategy end
function initial_condition(
    ::ZeroVelocityInitialCondition,
    q::Quadrature,
    problem::FluidFlowProblem,
    x::T,
    y::T
) where { T <: Real }
    copy(q.weights)
end
