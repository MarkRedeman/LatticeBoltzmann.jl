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
    q::Quadrature,
    x::Float64,
    y::Float64,
)
    nx = problem.NX
    ny = problem.NY
    nf = length(q.weights)

    f = Array{Float64}(undef, nx, ny, nf)

    x_range, y_range = range(problem)
    for x_idx = 1:nx, y_idx = 1:ny
        x = x_range[x_idx]
        y = y_range[y_idx]

        ρ = lattice_density(q, problem, x, y)
        u = lattice_velocity(q, problem, x, y),
        T = lattice_temperature(q, problem, x, y)

        f[x_idx, y_idx, :] = hermite_based_equilibrium(q, ρ, u, T)
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


initialize(
    ::InitializationStrategy,
    q::Quadrature,
    problem::FluidFlowProblem,
    cm = SRT
) = initialize(q, problem, cm)

function initialize(q::Quadrature, problem::FluidFlowProblem, cm = SRT)
    f = Array{Float64}(undef, problem.NX, problem.NY, length(q.weights))

    x_range, y_range = range(problem)
    for x_idx = 1:problem.NX, y_idx = 1:problem.NY
        f[x_idx, y_idx, :] =
            initial_condition(q, problem, x_range[x_idx], y_range[y_idx])
    end

    return f
end

function initial_condition(q::Quadrature, problem::FluidFlowProblem, x::Float64, y::Float64)

    ρ = lattice_density(q, problem, x, y)
    u = lattice_velocity(q, problem, x, y),
    T = lattice_temperature(q, problem, x, y)

    f = equilibrium(q, ρ, u, T)
    return f

    σ = deviatoric_tensor(q, problem, x, y, 0.0)
    for f_idx = 1:length(f)
        Q = hermite(Val{2}, q.abscissae[:, f_idx], q)
        f[f_idx] += q.weights[f_idx] * (cs^2 /factorial(2)) * dot(Q, σ)
    end
    return f

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
function initial_condition(q::Quadrature, problem::TGV, x::Float64, y::Float64)
    ρ = lattice_density(q, problem, x, y)
    u = lattice_velocity(q, problem, x, y)
    T = lattice_temperature(q, problem, x, y)
    f = equilibrium(q, ρ, u, T)

    τ = q.speed_of_sound_squared * lbm.lattice_viscosity(problem)
    # @show deviatoric_tensor(q, τ, f, ρ, u)

    cs = problem.q.speed_of_sound_squared
    τ = cs * problem.τ
    ν = problem.ν

    σ = deviatoric_tensor(q, problem, x, y, 0.0)
    ∇u = velocity_gradient(problem, x, y, 0.0)
    p = pressure(q, problem, x, y, 0.0)

    @show "HIER"
    @show σ ∇u (∇u + ∇u')
    τ = q.speed_of_sound_squared * lbm.lattice_viscosity(problem)

    @show σ + (∇u + ∇u') * viscosity(problem)


    @show viscosity(problem), ((τ + 0.5) * cs)

    for f_idx = 1:length(f)
        Q = hermite(Val{2}, q.abscissae[:, f_idx], q)
        # f[f_idx] += q.weights[f_idx] * τ * cs / (ν * factorial(2)) * dot(Q, σ)
        #
        # f[f_idx] += q.weights[f_idx] * (factor) * (cs /factorial(2)) * dot(Q, σ)


        # f[f_idx] += q.weights[f_idx] * (cs /factorial(2)) * dot(Q, (σ + δ * p) )

        # f[f_idx] += -q.weights[f_idx] * ((τ + 0.5) * cs) * dot(Q, ∇u)
        f[f_idx] += -q.weights[f_idx] * 0.5 * ((τ + 0.5) * cs) * dot(Q, ∇u + ∇u')
        # f[f_idx] += q.weights[f_idx] * dot(Q, σ) /2
        # @show Q
        # @show dot(Q, -∇u), dot(Q, σ)
        # @show dot(Q, -(∇u + ∇u') /2), dot(Q, σ)
        # f[f_idx] += q.weights[f_idx] * ((τ + 0.5) * cs) * dot(Q, σ) /2
    end

    # @show deviatoric_tensor(q, τ, f, ρ, u)
    # @show σ
    @show deviatoric_tensor(q, τ, f, ρ, u) ./ σ
    # @show σ
    return f
end
