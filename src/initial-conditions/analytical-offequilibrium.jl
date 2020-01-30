"""
1. Velocity, pressure and strain rate analytically
f₀(x) = f_eq(ρ₀(x), u₀(x), p₀(x)) + f¹(ρ₀(x), u₀(x), p₀(x), σ₀(x))
ρ₀ = ρ₀(x)
u₀ = u₀(x)
p₀ = p₀(x)
σ₀ = σ₀(x)
"""
struct AnalyticalEquilibriumAndOffEquilibrium <: InitializationStrategy end
function initial_condition(::AnalyticalEquilibriumAndOffEquilibrium, q::Quadrature, problem::FluidFlowProblem, x::T, y::T) where { T <: Real }
    f = equilibrium(q, problem, x, y)

    σ = problem.u_max^2 * deviatoric_tensor(q, problem, x, y, 0.0)
    ∇u = problem.u_max^2 * velocity_gradient(problem, x, y, 0.0)
    τ = q.speed_of_sound_squared * LatticeBoltzmann.lattice_viscosity(problem)

    cs = q.speed_of_sound_squared
    # @show σ ./ deviatoric_tensor(q, τ, f, density(q, f), velocity(q, f))

    # HMM..
    factor = (problem.domain_size[1] * problem.domain_size[2])

    # For Taylor Green Vortex...
    factor = 2 * (problem.domain_size[1] * problem.domain_size[2])

    for f_idx = 1:length(f)
        Q = hermite(Val{2}, q.abscissae[:, f_idx], q)

        # TODO: Check if we are using the correct viscosity here (we aren't...),
        # since we are setting the off equilibrium components of the shifted distribution funcitons
        f[f_idx] += - factor * q.weights[f_idx] * 0.5 * ((τ + 0.5) * cs) * dot(Q, ∇u + ∇u') #/ 0.785398163397454
        # f[f_idx] += q.weights[f_idx] * (cs^2 /factorial(2)) * dot(Q, σ)
    end

    velocity(q, f) = begin
        v = zeros(2)
        velocity!(q, f, density(q, f), v)
        v
    end
    σ = problem.u_max^2 * deviatoric_tensor(q, problem, x, y, 0.0)

    return f
end
function initial_condition(::AnalyticalEquilibriumAndOffEquilibrium, q::Quadrature, problem::TGV, x::T, y::T) where { T <: Real }
    @warn "Calling old TGV function"
    f = equilibrium(q, problem, x, y)

    σ = deviatoric_tensor(q, problem, x, y, 0.0)
    ∇u = velocity_gradient(problem, x, y, 0.0)
    p = pressure(q, problem, x, y, 0.0)

    @show "HIER"
    @show σ ∇u (∇u + ∇u')
    τ = q.speed_of_sound_squared * LatticeBoltzmann.lattice_viscosity(problem)
    cs = problem.q.speed_of_sound_squared
    τ = cs * problem.τ
    ν = problem.ν

    τ = q.speed_of_sound_squared * LatticeBoltzmann.lattice_viscosity(problem)

    @show σ + (∇u + ∇u') * viscosity(problem)
    @show viscosity(problem), ((τ + 0.5) * cs)

    for f_idx = 1:length(f)
        Q = hermite(Val{2}, q.abscissae[:, f_idx], q)
        # f[f_idx] += -q.weights[f_idx] * ((τ + 0.5) * cs) * dot(Q, ∇u)
        f[f_idx] += -q.weights[f_idx] * 0.5 * ((τ + 0.5) * cs) * dot(Q, ∇u + ∇u')
        # f[f_idx] += q.weights[f_idx] * dot(Q, σ) /2
        # f[f_idx] += q.weights[f_idx] * ((τ + 0.5) * cs) * dot(Q, σ) /2
    end

    return f
end
