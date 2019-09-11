using Test

@test 1 == 1

using lbm

quadratures = [
    D2Q4(),
    D2Q5(),
    D2Q9(),
    D2Q17(),
]
@testset "LBM.jl" begin

@testset "Multinomial hermite coefficients" begin

    ξ = [1.0 0.0]
    @test hermite(0, ξ) ≈ 1.0
    @test all(hermite(1, ξ) .≈ [1.0 0.0])
    @test all(hermite(2, ξ) .≈ [0.0 0.0; 0.0 -1.0])

    H_3 = hermite(3, ξ)
    @test all(H_3[:, :, 1] .≈ [-2.0  0.0;  0.0 -1.0])
    @test all(H_3[:, :, 2] .≈ [ 0.0 -1.0; -1.0  0.0])

    H_4 = hermite(4, ξ)

    # Check that our hermite polynomials obey the recurrance rule in 1 dimension
    ξ = 1.0
    for n = 3 : 4
        H = hermite(n, [ξ])[1]

        @test H ≈ ξ * hermite(n - 1, [ξ])[1] .- (n - 1) * hermite(n - 2, [ξ])[1]
    end

    # Check H4
    # H4 = [
    #     ξ[α_4] .* ξ[α_3] .* ξ[α_2] .* ξ[α_1] .- (
    #         # 6 terms
    #         (ξ[α_1] .* ξ[α_2]) .* δ(α_3, α_4) .+
    #         (ξ[α_1] .* ξ[α_3]) .* δ(α_2, α_4) .+
    #         (ξ[α_1] .* ξ[α_4]) .* δ(α_3, α_2) .+
    #         (ξ[α_2] .* ξ[α_3]) .* δ(α_1, α_4) .+
    #         (ξ[α_2] .* ξ[α_4]) .* δ(α_1, α_3) .+
    #         (ξ[α_3] .* ξ[α_4]) .* δ(α_1, α_2)
    #     ) .+ (δ(α_1, α_2) * δ(α_3, α_4) + δ(α_1, α_3) * δ(α_2, α_4) + δ(α_1, α_4) * δ(α_2, α_3))
    #     for α_1 = 1:D, α_2 = 1:D, α_3 = 1:D, α_4 = 1:D
    # ]
end

    include("quadrature.jl")

end
