abstract type CollisionModel end
struct SRT <: CollisionModel
    τ
end

function collide(cm::SRT, quadrature, f_in, τ)
    τ = collision_model.τ

    feq = equilibrium(quadrature, f_in);

    f_out = (1 - 1 / τ) * f_in + (1 / τ) * feq;

    return f_out
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
