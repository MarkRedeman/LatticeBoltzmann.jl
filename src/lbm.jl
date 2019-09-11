module lbm

include("quadratures/quadrature.jl")
include("stream.jl")
include("collision.jl")
# include("process.jl")
include("problems/problems.jl")

export Quadrature, Lattice, D2Q4, D2Q5, D2Q9, D2Q17
    initialize,
    density,
    momentum,
    total_energy,
    kinetic_energy,
    internal_energy,
    temperature,
    dimension,
    equilibrium,
    hermite_equilibrium,
    hermite_first_nonequilibrium

export stream
export CollisionModel,
    SRT,
    TRT,
    collide

# export process!

# abstract type Model

# struct Problem{Model}
#     model::Model
#     # model::Model
#     # quadrature::Quadrature
#     # collision::Collision
#     N::Int64      # rename to points / cells / ...
#     N_iter::Int64 # rename to iterations?
# end

# # Put collision_model, relaxation_rate, simulate, initial_condition, here..

# # abstract Quadrature;
# # Distribution{D2Q9}
# # Distribution{D2Q17}
# # immutable Distribution{T, Quadrature}
# #     fs::Array{T, 1}
# # end
# # typealias Distribution{T} Distribution{Float64, Q}

# typealias Distributions Array{Float64, 3}
# typealias Distribution Array{Float64, 2}

# include("quadratures/D2Q9.jl")
# # include("quadratures/D2Q4.jl")

# density(f::Distribution) = sum(f)

# velocity(f::Distribution) = velocity(f, abscissae, original_order)
# function velocity{D}(f::Distribution, abscissae::Array{Int64, D}, order::Array{Int64, 1})
#     ρ, u = density_and_velocity(f, abscissae, order)
#     return u
# end

# density_and_velocity(f::Distribution) = density_and_velocity(f, abscissae, original_order)
# function density_and_velocity{D}(f::Distribution, abscissae::Array{Int64, D}, order::Array{Int64, 1})::Tuple{Float64, Array{Float64, 1}}
#     u = zeros(D)
#     ρ = 0.0

#     # Compute: ∑fᵢ and ∑fᵢξᵢ
#     for idx ∈ order
#         for d = 1:D
#             u[d] += abscissae[d, idx] * f[idx]
#         end
#         ρ += f[idx]
#     end

#     return ρ, u / ρ
# end

# equilibrium(f::Distribution) = equilibrium(density_and_velocity(f, abscissae, original_order)...)
# equilibrium{T}(ρ::T, u::Array{T, 1})::Distribution = equilibrium(ρ, u, dot(u, u))
# equilibrium{T}(ρ::T, u::Array{T, 1}, u_squared::T)::Distribution = [equilibrium(ρ, u, u_squared, idx) for idx = 1:9]'
# function equilibrium{T}(rho::T, u::Array{T, 1}, u_squared::T, idx::Int)::T
#     const cs = dot(abscissae[:, idx], u)

#     return rho * weights[idx] .* (1.0 + 3.0 * cs + 4.5 * (cs .* cs) - 1.5 * u_squared)
# end

# """
# By default we will be using a bgk collision
# """
# collide(f::Distribution, ω)::Distribution = bgk_collision(f, ω)
# collide(f::Distribution, ω, force)::Distribution = bgk_collision(f, ω, force)

# """
# Apply the most simple collision operator
# """
# function bgk_collision{T}(f::Distribution, ω::T)::Distribution
#     const ρ, u = density_and_velocity(f)

#     bgk_collision(f, ω, ρ, u)
# end

# """
# Apply the bgk collision operator after adding additional momentum
# to the velocity due to a force term
# """
# function bgk_collision{T}(f::Distribution, ω::T, force::Array{T, 1})::Distribution
#     const ρ, u = density_and_velocity(f)

#     bgk_collision(f, ω, ρ, u + force / ω)
# end

# function bgk_collision{T}(f::Distribution, ω::T, ρ::T, u::Array{T, 1})
#     const u_squared = dot(u, u)

#     for idx = 1:9
#         f[idx] = (1 - ω) * f[ if next_x > lx
#         next_x -= lx
#     elseif next_x < 1
#         next_x += lx
#     end

#     next_y = y - abscissae[2, f_idx]
#     if next_y > ly
#         next_y -= ly
#     elseif next_y < 1
#         next_y += ly
#     end

#     return next_x, next_y
# end

δ(α, β) = α == β ? 1.0 : 0.0
hermite(::Type{Val{0}}) = 1.0
hermite(::Type{Val{0}}, ξ) = hermite(Val{0})
hermite(::Type{Val{1}}, ξ) = ξ .* hermite(Val{0})
function hermite(::Type{Val{2}}, ξ)
    D = length(ξ)
    return [
        ξ[β] * ξ[α] - δ(α, β) for α = 1:D, β = 1:D
    ]
end
function hermite(::Type{Val{3}}, ξ)
    D = length(ξ)
    return [
        ξ[γ] * ξ[β] * ξ[α] - (
            ξ[α] * δ(β, γ) + ξ[β] * δ(α, γ) + ξ[γ] * δ(α, β)
        )

        for α = 1:D, β = 1:D, γ = 1:D
    ]
end

function hermite(::Type{Val{4}}, ξ)
    N = 4
    D = length(ξ)
    H = fill(0.0, D, D, D, D)
    H_2 = hermite(Val{2}, ξ)
    H_3 = hermite(Val{3}, ξ)

    for i = 1 : D
        H[i, :, :, :] = ξ[i] * H_3

        for i_1 = 1:D, i_2 = 1:D, i_3 = 1:D
            H[i, i_1, i_2, i_3] -= δ(i, i_1) * H_2[i_2, i_3] - δ(i, i_2) * H_2[i_1, i_3] - δ(i, i_3) * H_2[i_1, i_2]
        end
    end
    return H
end

# hermite(::Type{Val{N}}, ξ) where N = begin
#     D = length(ξ)
#     H = fill(0.0, [D for 1:N]...)

#     H_n_1 = hermite(Val{N - 1}, ξ)
#     H_n_2 = hermite(Val{N - 2}, ξ)
#     for i = 1 : D
#         H_d = sum([
#             δ(i, i_k) *
#         ])
#         H[i, :] = ξ[i] * H_n_1 H_d
#     end

#     hermite(Val{N - 1}, ξ)
# end
hermite(N::Int, ξ) = hermite(Val{N}, ξ)

using TensorOperations
function hermite(::Type{Val{5}}, x)
    D = length(x)
    N = 5
    H = fill(undef, [D for _ in 1:N]...)
    H_4 = hermite(Val{4}, x)
    H_3 = hermite(Val{3}, x)

    δ(α, β) = α == β ? 1.0 : 0.0

    @tensor begin
        H[i, i1, i2, i3, i4] = x[i] * H[i1, i2, i3, i4] - (
            δ(i, i1) * H_3[i2, i3, i4] +
            δ(i, i2) * H_3[i1, i3, i4] +
            δ(i, i3) * H_3[i1, i2, i4] +
            δ(i, i4) * H_3[i1, i2, i3]
        )
    end
end

export hermite

end
