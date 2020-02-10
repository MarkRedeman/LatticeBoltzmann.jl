"""
TODO
"""
struct AnalyticalVelocityAndStress <: InitializationStrategy end
function initial_condition(::AnalyticalVelocityAndStress, q::Quadrature, problem::FluidFlowProblem, x::T, y::T) where { T <: Real }
    u = lattice_velocity(q, problem, x, y)
    ρ = 1.0

    f = equilibrium(q, one(T), u, one(T))

    σ = problem.u_max^2 * deviatoric_tensor(q, problem, x, y, 0.0)
    ∇u = problem.u_max^2 * velocity_gradient(problem, x, y, 0.0)

    # for f_idx = 1:length(f)
    #     Q = hermite(Val{2}, q.abscissae[:, f_idx], q)

    #     # TODO: Check if we are using the correct viscosity here (we aren't...),
    #     # since we are setting the off equilibrium components of the shifted distribution funcitons
    #     f[f_idx] += - factor * q.weights[f_idx] * 0.5 * ((τ + 0.5) * cs) * dot(Q, ∇u + ∇u') #/ 0.785398163397454
    #     # f[f_idx] += q.weights[f_idx] * (cs^2 /factorial(2)) * dot(Q, σ)
    # end

    cs = q.speed_of_sound_squared
    τ = cs * lattice_viscosity(problem)
    τ_eff = τ + 0.5
    for f_idx = 1:length(f)
        w_i = q.weights[f_idx]
        H_2 = hermite(Val{2}, q.abscissae[:, f_idx], q)


        # TODO: possibly replace with pressure
        temperature = 1.0
        f_neq = w_i * (cs * τ_eff * ρ * temperature) / (2) * dot(H_2, ∇u + ∇u')
        f[f_idx] += - f_neq
    end


    @info "Check"
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
end
