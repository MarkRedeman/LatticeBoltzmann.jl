abstract type CollisionModel end
struct SRT <: CollisionModel
    τ
end

function collide(collision_model::SRT, q, f_in)
    # Density
    ρ = density(q, f_in)

    # Momentum
    j = momentum(q, f_in)

    # Temperature
    T = 1.0
    # T = temperature(q, f_in, ρ, j ./ ρ)

    feq = equilibrium(q, ρ, j ./ ρ, T);

    τ = collision_model.τ
    f_out = (1 - 1 / τ) * f_in + (1 / τ) * feq;


    # Since this is the place where we will have computed all
    # hermite coefficients, call process! here? (or at least a callback)
    # cb(f_out, moments)

    return f_out
end

struct SRT_Force <: CollisionModel
    τ
    F
end

struct TRT <: CollisionModel
    τ_symmetric
    τ_asymmetric
end

function collide(collision_model::TRT, quadrature, f_in)
    τ_s = collision_model.τ_symmetric
    τ_a = collision_model.τ_asymmetric

    # TODO
    return f_in
end