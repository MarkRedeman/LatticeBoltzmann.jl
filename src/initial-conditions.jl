abstract type InitializationStrategy end
InitializationStrategy(problem) = FromEquilibrium2()

# 1. Velocity, pressure and strain rate analytically
# 2. Velocity and pressure only (initialize with equilibrium)
# 3. Veloctiy and stress only (ρ_0 = 1.0)
# 4. Initialization scheme from Mei et al

# IDEA: use traits to determine if a problem is a steady state problem
# if so, we can return a zero initial condition,
# otherwise we want to determine the condition based on availability
# of derivatives

struct ZeroInitialCondition2 <: InitializationStrategy end

struct FromEquilibrium2 <: InitializationStrategy end

struct WithNonEquilibrium2 <: InitializationStrategy end

struct MeiEtAl2 <: InitializationStrategy end

function initialize(
    strategy::InitializationStrategy,
    problem::FluidFlowProblem,
    quadrature::Quadrature,
    x::Float64,
    y::Float64,
)
    nx = problem.NX
    ny = problem.NY
    nf = length(quadrature.weights)

    f = Array{Float64}(undef, nx, ny, nf)

    x_range, y_range = range(problem)
    for x_idx = 1:nx, y_idx = 1:ny
        x = x_range[x_idx]
        y = y_range[y_idx]

        ρ =
            lattice_density(quadrature, problem, x, y), u =
                lattice_velocity(quadrature, problem, x, y), T =
                    lattice_temperature(quadrature, problem, x, y)

        f[x_idx, y_idx, :] = hermite_based_equilibrium(quadrature, ρ, u, T)
    end

    return f
end

# function initial_equilibrium(quadrature::Quadrature, problem::FluidFlowProblem, x::Float64, y::Float64)
#     return equilibrium(
#         quadrature,
#         lattice_density(quadrature, problem, x, y),
#         lattice_velocity(quadrature, problem, x, y),
#         lattice_temperature(quadrature, problem, x, y)
#     )
# end

# function initial_condition(q::Quadrature, problem::FluidFlowProblem, x::Float64, y::Float64)
#     initial_equilibrium(q, problem, x, y)
# end

# function initialize(quadrature::Quadrature, problem::FluidFlowProblem, cm = SRT)
#     f = Array{Float64}(undef, problem.NX, problem.NY, length(quadrature.weights))

#     x_range, y_range = range(problem)
#     for x_idx in 1:problem.NX, y_idx in 1:problem.NY
#         f[x_idx, y_idx, :] = initial_equilibrium(
#             quadrature,
#             problem,
#             x_range[x_idx],
#             y_range[y_idx]
#         )
#     end

#     return f
# end
