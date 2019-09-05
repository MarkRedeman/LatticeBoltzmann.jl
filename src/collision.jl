abstract type CollisionModel end
struct SRT <: CollisionModel
    τ
end

function collide(collision_model::SRT, quadrature, f_in)
    τ = collision_model.τ

    feq = equilibrium(quadrature, f_in);

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

# """
# Apply the most simple collision operator
# """
# function bgk_collision{T}(f::Distribution, ω::T)::Distribution
#     ρ, u = density_and_velocity(f)

#     bgk_collision(f, ω, ρ, u)
# end

# """
# Apply the bgk collision operator after adding additional momentum
# to the velocity due to a force term
# """
# function bgk_collision{T}(f::Distribution, ω::T, force::Array{T, 1})::Distribution
#     ρ, u = density_and_velocity(f)

#     bgk_collision(f, ω, ρ, u + force / ω)
# end

# function bgk_collision{T}(f::Distribution, ω::T, ρ::T, u::Array{T, 1})
#     u_squared = dot(u, u)

#     for idx = 1:9
#         f[idx] = (1 - ω) * f[idx] + ω * equilibrium(ρ, u, u_squared, idx)
#     end

#     return f
# end
