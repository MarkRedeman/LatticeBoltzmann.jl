abstract type InitializationStrategy end

# TODO: base this on the existence of a velocity gradient, pressure and velocity
InitializationStrategy(problem) = AnalyticalEquilibriumAndOffEquilibrium()

function initialize(strategy::InitializationStrategy, q::Quadrature, problem::FluidFlowProblem, cm = SRT)
    f = Array{Float64}(undef, problem.NX, problem.NY, length(q.weights))

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
    x,
    y
) = initial_condition(q, problem, x, y)


# 1. Velocity, pressure and strain rate analytically
"""
f₀(x) = f_eq(ρ₀(x), u₀(x), p₀(x)) + f¹(ρ₀(x), u₀(x), p₀(x), σ₀(x))
ρ₀ = ρ₀(x)
u₀ = u₀(x)
p₀ = p₀(x)
σ₀ = σ₀(x)
"""
struct AnalyticalEquilibriumAndOffEquilibrium <: InitializationStrategy end
function initial_condition(::AnalyticalEquilibriumAndOffEquilibrium, q::Quadrature, problem::FluidFlowProblem, x::Float64, y::Float64)
    ρ = lattice_density(q, problem, x, y)
    u = lattice_velocity(q, problem, x, y)
    T = lattice_temperature(q, problem, x, y)
    f = equilibrium(q, ρ, u, T)

    σ = deviatoric_tensor(q, problem, x, y, 0.0)

    cs = q.speed_of_sound_squared

    for f_idx = 1:length(f)
        Q = hermite(Val{2}, q.abscissae[:, f_idx], q)
        f[f_idx] += q.weights[f_idx] * (cs^2 /factorial(2)) * dot(Q, σ)
    end
    return f
end
function initial_condition(::AnalyticalEquilibriumAndOffEquilibrium, q::Quadrature, problem::TGV, x::Float64, y::Float64)
    @warn "Calling old TGV function"
    ρ = lattice_density(q, problem, x, y)
    u = lattice_velocity(q, problem, x, y)
    T = lattice_temperature(q, problem, x, y)
    f = equilibrium(q, ρ, u, T)

    σ = deviatoric_tensor(q, problem, x, y, 0.0)
    ∇u = velocity_gradient(problem, x, y, 0.0)
    p = pressure(q, problem, x, y, 0.0)

    @show "HIER"
    @show σ ∇u (∇u + ∇u')
    τ = q.speed_of_sound_squared * lbm.lattice_viscosity(problem)
    cs = problem.q.speed_of_sound_squared
    τ = cs * problem.τ
    ν = problem.ν

    τ = q.speed_of_sound_squared * lbm.lattice_viscosity(problem)

    @show σ + (∇u + ∇u') * viscosity(problem)
    @show viscosity(problem), ((τ + 0.5) * cs)

    for f_idx = 1:length(f)
        Q = hermite(Val{2}, q.abscissae[:, f_idx], q)
        # f[f_idx] += -q.weights[f_idx] * ((τ + 0.5) * cs) * dot(Q, ∇u)
        f[f_idx] += -q.weights[f_idx] * 0.5 * ((τ + 0.5) * cs) * dot(Q, ∇u + ∇u')
        # f[f_idx] += q.weights[f_idx] * dot(Q, σ) /2
        # f[f_idx] += q.weights[f_idx] * ((τ + 0.5) * cs) * dot(Q, σ) /2
    end

    # @show deviatoric_tensor(q, τ, f, ρ, u)
    # @show σ
    @show deviatoric_tensor(q, τ, f, ρ, u) ./ σ
    return f
end

# 2. Velocity and pressure only (initialize with equilibrium)
"""
f₀(x) = f_eq(ρ₀(x), u₀(x), p₀(x))
ρ₀ = ρ₀(x)
u₀ = u₀(x)
p₀ = p₀(x)
"""
struct AnalyticalEquilibrium <: InitializationStrategy end
function initial_condition(::AnalyticalEquilibrium, q::Quadrature, problem::FluidFlowProblem, x::Float64, y::Float64)
    ρ = lattice_density(q, problem, x, y)
    u = lattice_velocity(q, problem, x, y)
    T = lattice_temperature(q, problem, x, y)

    return equilibrium(q, ρ, u, T)
end

# 3. Veloctiy and stress only (ρ_0 = 1.0)
"""
f₀(x) = f_eq(ρ₀(x), u₀(x), p₀(x))
ρ₀ = 1.0
u₀ = u₀(x)
p₀ = p₀(x)
"""
struct AnalyticalVelocity <: InitializationStrategy end

# TODO compute the density
average_density(problem::FluidFlowProblem) = 1.0
average_pressure(problem::FluidFlowProblem) = 1.0
function density(
    q::Quadrature,
    problem::FluidFlowProblem,
    x,
    y,
    p̄ = average_pressure(problem),
    ρ̄ = average_density(problem)
)
    δρ = (pressure(problem, x, y) - p̄) / q.speed_of_sound_squared

    return ρ̄ + δρ
end
# TODO Move computing the averages to the initialize function so that these can be computed globaly
function initial_condition(::AnalyticalVelocity, q::Quadrature, problem::FluidFlowProblem, x::Float64, y::Float64)
    p̄ = average_pressure(problem)
    ρ̄ = average_density(problem)

    # Compute the density from the pressure, average pressure and average density
    ρ = density(q, problem, x, y, p̄, ρ̄)

    u = lattice_velocity(q, problem, x, y),
    T = pressure(problem, x, y) / ρ

    return equilibrium(q, ρ, u, T)
end

# 4. Initialization scheme from Mei et al
struct IterativeInitialization <: InitializationStrategy
    τ
end
IterativeInitialization() = IterativeInitialization(0.8)
