
struct D2Q9 <: Quadrature end
# type D2Q9 <: Lattice end

# weights(lattice::D2Q9) = [4/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9]'
# abscissaes(lattice::D2Q9) = [0  0; -1  1; -1  0; -1 -1; 0 -1; 1 -1; 1  0; 1  1; 0  1]'
# opposite(lattice::D2Q9) = [1 2 3 4 5 6 7 8 9]
# opposite(lattice::D2Q9, idx) = idx + (9 - 1) / 2
# dimension(lattice::D2Q9) = 2
# order(lattice::D2Q9) = [1, 7, 9, 3, 5, 8, 2, 4, 6]

# D2Q9
const original_order = [1, 7, 9, 3, 5, 8, 2, 4, 6]
const weights   = [4/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9]'
const abscissae = [
    0  -1  -1  -1   0   1  1  1  0
    0   1   0  -1  -1  -1  0  1  1
]


const speed_of_sound = 1 / √3

struct D2Q9 <: Quadrature end

density(::D2Q9, f) = sum(f, dims=3)
function momentum(::D2Q9, f)
    jx = sum(f[:,:,[6, 7, 8]], dims=3) .- sum(f[:,:,[2, 3, 4]], dims=3);
    jy = sum(f[:,:,[2, 8, 9]], dims=3) .- sum(f[:,:,[4, 5, 6]], dims=3);

    jx, jy
end

function equilibrium(quadrature::D2Q9, f_in)
    # Density
    f_ρ = density(quadrature, f_in)

    # Momentum
    jx, jy = momentum(quadrature, f_in)

    # Velocity componetns
    u = jx ./ f_ρ
    v = jy ./ f_ρ

    feq = equilibrium(f_ρ, u, v);

    @show sum(f_in), sum(feq), sum(jx + jy), sum(jx.^2 .+ jy.^2)
    return feq
    # @show u.^2 .+ v.^2
    # @show u

    # Compute equilibrium distribution
end

function equilibrium!(::D2Q9, f, ρ, velocity, T)
    u = velocity[1]
    v = velocity[2]

    u_squared = u.^2 .+ v.^2

    for idx = 1:9
        cs = abscissae[1, idx] .* u .+ abscissae[2, idx] .* v


        f[idx] = ρ .* weights[idx] .* (1.0 .+ 3.0 * cs .+ 4.5 * (cs .* cs) .- 1.5 * u_squared)
    end

    return f
end


function equilibrium(q::D2Q9, ρ, velocity, T)
    f = fill(0.0, 9)
    equilibrium!(q, f, ρ, velocity, T)
    return f
end


function equilibrium(rho, u, u_squared, idx::Int)
    @info "HUH"
    cs = dot(abscissae[:, idx], u)

    return rho * weights[idx] .* (1.0 + 3.0 * cs + 4.5 * (cs .* cs) - 1.5 * u_squared)
end

function equilibrium(ρ,u,v)
    f = zeros(size(ρ,1), size(ρ,2), 9);

    # subexp1 = 1 .- 1.5*(u.^2+v.^2);

    # f[:, :, 1] = (1/9) .*  ρ .* (subexp1 + 3 .* u .+ 9 ./ 2 .* u.^2);
    # f[:, :, 2] = (1/9) .*  ρ .* (subexp1 + 3 .* v .+ 9 ./ 2 .* v.^2);
    # f[:, :, 3] = (1/9) .*  ρ .* (subexp1 - 3 .* u .+ 9 ./ 2 .* u.^2);
    # f[:, :, 4] = (1/9) .*  ρ .* (subexp1 - 3 .* v .+ 9 ./ 2 .* v.^2);

    # f[:, :, 5] = (1/36) .* ρ .* (subexp1 + 3 .* ( u .+ v) .+ 9 ./ 2 .* ( u .+ v).^2);
    # f[:, :, 6] = (1/36) .* ρ .* (subexp1 + 3 .* (-u .+ v) .+ 9 ./ 2 .* (-u .+ v).^2);
    # f[:, :, 7] = (1/36) .* ρ .* (subexp1 + 3 .* (-u .- v) .+ 9 ./ 2 .* (-u .- v).^2);
    # f[:, :, 8] = (1/36) .* ρ .* (subexp1 + 3 .* ( u .- v) .+ 9 ./ 2 .* ( u .- v).^2);

    # f[:, :, 9] = (4/9) .* ρ .* (subexp1);

    u_squared = u.^2 + v.^2
    for idx = 1:9
        cs = abscissae[1, idx] .* u .+ abscissae[2, idx] .* v


        f[:, :, idx] = ρ .* weights[idx] .* (1.0 .+ 3.0 * cs .+ 4.5 * (cs .* cs) .- 1.5 * u_squared)
    end

    return f
end


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

# equilibrium(set::D2Q9, f) = 1.
# equilibrium(set::D2Q9, ρ, u) = 1.
# equilibrium(set::D2Q9, ρ, u, θ) = 1.

# function equilibrium{T}(rho::T, u::Array{T, 1}, u_squared, idx)::T
#     const cs = dot(abscissae[:, idx], u)

#     return rho * weights[idx] .* (1.0 + 3.0 * cs + 4.5 * (cs .* cs) - 1.5 * u_squared)
# end

# Note: check "Galilaen invariance of lbm" for a nice overview
# Q
# quadrature_precission(::Type{D2Q9}) = 5
# weights(::Type{D2Q9}) = [4/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9]'
# abscissae(::Type{D2Q9}) = [
#     0  -1  -1  -1   0   1  1  1  0
#     0   1   0  -1  -1  -1  0  1  1
# ]
# speed_of_sound(::Type{D2Q9}) = 1 / √3
# equilibrium{N<:Int} = equilibrium{2}
# equilibrium{1} = 1
# equilibrium{2} = 2


# abstract type Lattice end

# struct D2Q9 <: Lattice
#     const weights = [4/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9]
#     const velocities = [
#         0  -1  -1  -1   0   1  1  1  0
#         0   1   0  -1  -1  -1  0  1  1
#     ]
#     const Q = length(weights)
#     const opposites = [1 6 7 8 9 2 3 4 5]
# end
