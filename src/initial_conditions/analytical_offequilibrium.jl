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
    factor = (problem.domain_size[1] * problem.domain_size[2])

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
    # @warn "Calling old TGV function"
    ρ = lattice_density(q, problem, x, y)
    f = equilibrium(
        q,
        ρ,
        lattice_velocity(q, problem, x, y),
        lattice_temperature(q, problem, x, y)
    )

    u = lattice_velocity(q, problem, x, y)
    # ρ = 1.0
    # f = equilibrium(q, one(T), u, one(T))

    # ρ = 1.0
    σ = problem.u_max^2 * deviatoric_tensor(q, problem, x, y, 0.0)
    ∇u = problem.u_max^2 * velocity_gradient(problem, x, y, 0.0)


    σ = deviatoric_tensor(q, problem, x, y, 0.0)
    ∇u = velocity_gradient(problem, x, y, 0.0)
    p = pressure(q, problem, x, y, 0.0)

    # NEW
    cs = q.speed_of_sound_squared
    τ = cs * lattice_viscosity(problem)
    τ_eff = τ + 0.5
    for f_idx = 1:length(f)
        w_i = q.weights[f_idx]
        H_2 = hermite(Val{2}, q.abscissae[:, f_idx], q)


        # TODO: possibly replace with pressure
        temperature = 1.0
        f_neq = w_i * (cs * τ_eff * ρ * temperature) / (2) * dot(H_2, ∇u + ∇u')

        # ν = τ_eff / cs
        # f_neq = w_i * (cs^2 * ν) / (2) * dot(H_2, ∇u + ∇u')

        f[f_idx] -= f_neq
    end


    # @info "Check"
    u = lattice_velocity(q, problem, x, y)
    u_lb = copy(u)
    ρ = density(q, f)
    velocity!(q, f, ρ, u_lb)
    # @show u ./ u_lb

    σ = deviatoric_tensor(q, problem, x, y, 0.0)
    σ_lb = deviatoric_tensor(q, τ, f, ρ, u)
    # @show σ ./ σ_lb
    ρ_lb = ρ
    ρ = lattice_density(q, problem, x, y)
    # @show ρ ./ ρ_lb





    return f


    @show "HIER"
    # @show σ ∇u (∇u + ∇u')
    τ = q.speed_of_sound_squared * LatticeBoltzmann.lattice_viscosity(problem)
    cs = problem.q.speed_of_sound_squared
    τ = cs * problem.τ
    ν = problem.ν

    τ = q.speed_of_sound_squared * LatticeBoltzmann.lattice_viscosity(problem)

    # @show σ + (∇u + ∇u') * viscosity(problem)
    # @show viscosity(problem), ((τ + 0.5) * cs)


    # TODO: Check if density should be taken into account...
    # Λ = ();

    @show "Testing"
    for f_idx = 1:length(f)
        Q = hermite(Val{2}, q.abscissae[:, f_idx], q)
        # f[f_idx] += -q.weights[f_idx] * ((τ + 0.5) * cs) * dot(Q, ∇u)
        f[f_idx] += - ρ * q.weights[f_idx] * 0.5 * ((τ + 0.5) * cs) * dot(Q, ∇u + ∇u')
        # f[f_idx] += q.weights[f_idx] * dot(Q, σ) /2
        # f[f_idx] += q.weights[f_idx] * ((τ + 0.5) * cs) * dot(Q, σ) /2
    end

    return f
end
