# functions to determine hermite coefficients of a function

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
