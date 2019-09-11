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

@testset "Symmetry requirements of lattice $q" for q in quadratures
    w = q.weights
    ξs = q.abscissae
    D = 2
    Q = length(w)
    δ(α, β) = α == β ? 1.0 : 0.0

    @testset "Symmetry $n" for n = 0:lbm.order(q)
        if n == 0
            H = sum(w) - 1.0
        elseif n == 1
            H = [
                sum(ξs[α, :] .* w[:])
                for α = 1:D
            ]
        elseif n == 2
            H = [
                sum(w[:] .* (q.speed_of_sound_squared .* ξs[β, :] .* ξs[α, :] .- δ(α, β)))
                for α = 1:D, β = 1:D
            ]
        elseif n == 3
            H = [
                sum(
                    w[:] .* (
                        (q.speed_of_sound_squared .* ξs[γ, :] .* ξs[β, :] .* ξs[α, :]) .-
                        (ξs[α, :] .* δ(β, γ) .+ ξs[β, :] .* δ(α, γ) .+ ξs[γ, :] .* δ(α, β))
                    )
                )
                for α = 1:D, β = 1:D, γ = 1:D
            ]
        elseif n == 4
            H = [
                sum(
                    w[:] .* (
                        q.speed_of_sound_squared^2 .* ξs[α_4, :] .* ξs[α_3, :] .* ξs[α_2, :] .* ξs[α_1, :] .-
                        q.speed_of_sound_squared .* (
                            # 6 terms
                            (ξs[α_1, :] .* ξs[α_2, :]) .* δ(α_3, α_4) .+
                            (ξs[α_1, :] .* ξs[α_3, :]) .* δ(α_2, α_4) .+
                            (ξs[α_1, :] .* ξs[α_4, :]) .* δ(α_3, α_2) .+
                            (ξs[α_2, :] .* ξs[α_3, :]) .* δ(α_1, α_4) .+
                            (ξs[α_2, :] .* ξs[α_4, :]) .* δ(α_1, α_3) .+
                            (ξs[α_3, :] .* ξs[α_4, :]) .* δ(α_1, α_2)
                        ) .+
                        (δ(α_1, α_2) * δ(α_3, α_4) + δ(α_1, α_3) * δ(α_2, α_4) + δ(α_1, α_4) * δ(α_2, α_3))
                    )
                )
                for α_1 = 1:D, α_2 = 1:D, α_3 = 1:D, α_4 = 1:D
            ]
        elseif n == 5
            H = 0.0
            # H = [
            #     sum(
            #         w[:] .* (
            #             q.speed_of_sound_squared^2 .* ξs[α_5, :] .* ξs[α_4, :] .* ξs[α_3, :] .* ξs[α_2, :] .* ξs[α_1, :] .-
            #             q.speed_of_sound_squared .* (
            #                 # 10 terms (5!) / (2^2 * 3! * 1!)
            #                 (ξs[α_1, :] .* ξs[α_2, :] .* ξs[α_3, :]) .* δ(α_4, α_5) .+
            #             ) .+
            #             (
            #                 # 5! / (2^2 1! 2!) = 15 terms :O
            #                 δ(α_1, α_2) * δ(α_3, α_4) .* ξ[α_5, :] .+
            #             )
            #         )
            #     )
            #     for α_1 = 1:D, α_2 = 1:D, α_3 = 1:D, α_4 = 1:D, α_5 = 1:D
            # ]
        else
            H = 0.0
        end

        @test all(abs.(H) .< 10eps())
    end
end

@testset "moments of the equilibrium distribution of the $q lattice" for q in quadratures
end

end
