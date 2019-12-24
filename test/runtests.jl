using Test

@test 1 == 1

using lbm

quadratures = [D2Q4(), D2Q5(), D2Q9(), D2Q13(), D2Q17(), D2Q21(), D2Q37()]
quadratures = lbm.Quadratures

include("problems/taylor-green-vortex.jl")
include("problems/poiseuille.jl")

@testset "LBM.jl" begin

    @testset "Multinomial hermite coefficients" begin

        ξ = [1.0 0.0]
        @test hermite(0, ξ) ≈ 1.0
        @test all(hermite(1, ξ) .≈ [1.0 0.0])
        @test all(hermite(2, ξ) .≈ [0.0 0.0; 0.0 -1.0])

        H_3 = hermite(3, ξ)
        @test all(H_3[:, :, 1] .≈ [-2.0 0.0; 0.0 -1.0])
        @test all(H_3[:, :, 2] .≈ [0.0 -1.0; -1.0 0.0])

        H_4 = hermite(4, ξ)

    # Check that our hermite polynomials obey the recurrance rule in 1 dimension
        ξ = 1.0
        for n = 3:4
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



    logspace(start, stop, length) = exp10.(range(start, stop = stop, length = length))
    @testset "moments of the equilibrium distribution of the $q lattice" for q in quadratures
        N = 1

        f = lbm.equilibrium(q, 1.0, zeros(lbm.dimension(q)), 1.0)

        @test lbm.density(q, f) ≈ 1.0

        v = zeros(lbm.dimension(q))
        lbm.velocity!(q, f, 1.0, v)
        @test isapprox(v, zeros(lbm.dimension(q)), atol = 1e-16)

        equilibrium(ρ, u) = begin
            T = 1.0
            # return lbm.equilibrium(q, ρ, u, T)

            f = zeros(length(q.weights))
            lbm.hermite_based_equilibrium!(q, ρ, u, T, f)
            f1 = copy(f)
            lbm.equilibrium!(q, ρ, u, T, f)

            if ! isapprox(lbm.density(q, f), 1.0, rtol = 1e-10, atol = 1e-8)
                # @show f1 f
                # @show f1 - f
                lbm.hermite_based_equilibrium!(q, ρ, u, T, f)
            end
            if ! isapprox(lbm.density(q, f1), 1.0, rtol = 1e-10, atol = 1e-8)
                # @show f1 f
                # @show f1 - f
                lbm.hermite_based_equilibrium!(q, ρ, u, T, f)
            end
            # @show isapprox(lbm.density(q, f), 1.0, rtol = 1e-10, atol = 1e-8) isapprox(lbm.density(q, f1), 1.0, rtol = 1e-10, atol = 1e-8)
            # @show f - f1
            # return f1
            return f
        end
        density(ρ, u) = lbm.density(q, equilibrium(ρ, u))
        momentum(ρ, u) = begin
            v = zeros(lbm.dimension(q))
            lbm.velocity!(q, equilibrium(ρ, u), 1.0, v)
            return v
        end

        pressure(ρ, u) = LBM.pressure(equilibrium(ρ, u))

        @test density(1.0, [0.0, 0.0]) ≈ 1.0
        @test isapprox(
            momentum(1.0, [0.0, 0.0]),
            fill(0.0, lbm.dimension(q)),
            atol = 1e-16,
        )

    # Todo: analytically check the numerical error bounds

        for u ∈ logspace(2, -9, 12)
            for v ∈ [[0, 0], [u, 0.0], [0.0, u], [u, u]]
                @test isapprox(density(1.0, v), 1.0, rtol = 1e-10, atol = 1e-8)
                for d = 1:lbm.dimension(q)
                    @test isapprox(
                        momentum(1.0, v)[d],
                        v[d],
                        atol = 1e-11,
                        rtol = 1e-6,
                    )
                end
            end
        end
    end
end
