struct IterativeInitializationCollisionModel{T <: Real, UT <: AbstractArray{T, 3}} <:
       CollisionModel
    τ::T

    # The incompressible equilibrium funciton can be split in a linear and a nonlinear component
    # the linear component depends on the local density while the nonlinear component depends
    # on the local velocity.
    # Since this collision operator keeps a constant velocity for each lattice node we can
    # precompute the nonlinear term for some performance improvements
    nonlinear_term::UT
end
function IterativeInitializationCollisionModel(
    q::Quadrature,
    τ::T,
    problem,
) where {T <: Real}
    nx, ny, nf = problem.NX, problem.NY, length(q.weights)
    nonlinear_term = Array{T}(undef, nx, ny, nf)
    ρ_0 = one(T)

    # Initialize f with the given velocity field
    x_range, y_range = range(problem)
    for x_idx in 1:nx, y_idx in 1:ny
        u = lattice_velocity(q, problem, x_range[x_idx], y_range[y_idx])
        u_squared = u[1]^2 + u[2]^2

        for f_idx in 1:nf
            u_dot_xi = dot(q.abscissae[:, f_idx], u)

            cs = q.speed_of_sound_squared
            a_H_1 = cs * u_dot_xi
            a_H_2 = cs^2 * (u_dot_xi * u_dot_xi) - cs * u_squared

            nonlinear_term[x_idx, y_idx, f_idx] =
                q.weights[f_idx] * ρ_0 * (a_H_1 + a_H_2 / 2)
        end
    end

    return IterativeInitializationCollisionModel(τ, nonlinear_term)
end
function collide!(
    collision_model::IterativeInitializationCollisionModel,
    q::Quadrature,
    f_in,
    f_out;
    time = 0.0,
)
    τ = collision_model.τ
    nx, ny, nf = size(f_out)
    for x_idx in 1:nx, y_idx in 1:ny
        ρ = density(q, f_in[x_idx, y_idx, :])
        @inbounds for f_idx in 1:nf
            feq = q.weights[f_idx] * ρ + collision_model.nonlinear_term[x_idx, y_idx, f_idx]

            f_out[x_idx, y_idx, f_idx] =
                (1 - 1 / τ) * f_in[x_idx, y_idx, f_idx] + (1 / τ) * feq
        end
    end
    return
end
