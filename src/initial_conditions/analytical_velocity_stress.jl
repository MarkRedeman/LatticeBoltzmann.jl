"""
TODO
"""
struct AnalyticalVelocityAndStress <: InitializationStrategy end
function initial_condition(
    ::AnalyticalVelocityAndStress,
    q::Quadrature,
    problem::FluidFlowProblem,
    x::T,
    y::T,
) where {T <: Real}
    u = lattice_velocity(q, problem, x, y)
    ρ = 1.0
    f = equilibrium(q, one(T), u, one(T))
    σ = problem.u_max^2 * deviatoric_tensor(q, problem, x, y, 0.0)
    ∇u = problem.u_max^2 * velocity_gradient(problem, x, y, 0.0)

    cs = q.speed_of_sound_squared
    τ = cs * lattice_viscosity(problem)
    τ_eff = τ + 0.5
    for f_idx in 1:length(f)
        w_i = q.weights[f_idx]
        H_2 = hermite(Val{2}, q.abscissae[:, f_idx], q)

        # TODO: possibly replace with pressure
        temperature = 1.0
        f_neq = w_i * (cs * τ_eff * ρ * temperature) / (2) * dot(H_2, ∇u + ∇u')
        f[f_idx] += -f_neq
    end
    return f
end
