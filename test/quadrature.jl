using BenchmarkTools
using TimerOutputs

@testset "Symmetry requirements of lattice $q" for q in quadratures
    w = q.weights
    ξs = q.abscissae
    D = 2
    Q = length(w)
    δ(α, β) = α == β ? 1.0 : 0.0

    # If we didn't manage to choose correct weights all other tests will likely fail
    @test sum(q.weights) ≈ 1.0

    @testset "Symmetry $n" for n = 1:div(lbm.order(q), 2)
        H = sum([
            q.weights[f_idx] *
            hermite(n, q.abscissae[:, f_idx], q)
            for f_idx = 1:length(q.weights)
        ])

        @test all(isapprox.(H, 0.0, atol=1e-15))
    end
end

@testset "moments of the equilibrium distribution of the $q lattice" for q in quadratures

    N = 1
    ρ = fill(1.0, N, N)
    u = fill(0.001, N, N, lbm.dimension(q))
    T = fill(1.0, N, N)

    @testset "equilibira" begin
        f = lbm.equilibrium(q, ρ, u, T);

        @test all(lbm.density(q, f) .≈ ρ)
        @test all(lbm.momentum(q, f) ./ ρ .≈ u)
        @test all(isapprox.(lbm.temperature(q, f, ρ, u), T, atol=1e-5))
        # @code_warntype lbm.equilibrium(q, ρ, u, T)
        @inferred lbm.equilibrium(q, ρ, u, T)

    end

    @testset "temperature" begin
        N = 1
        ρ = fill(1.0, N, N)
        u = fill(0.0, N, N, lbm.dimension(q))
        T = fill(1.0, N, N)
        f = lbm.equilibrium(q, ρ, u, T);
        # @show lbm.temperature(q, f, ρ, u) * q.speed_of_sound_squared
        @test all(lbm.temperature(q, f, ρ, u) .≈ T)
    end

    @testset "moments of equilibira in a single point" begin
        ρ = 1.0
        u = [0.1, 0.1]
        T = 1.0
        f = lbm.equilibrium(q, ρ, u, T);

        @test lbm.density(q, f) ≈ ρ
        @test all(lbm.momentum(q, f) ./ ρ .≈ u)

        @inferred lbm.density(q, f)
        @inferred lbm.momentum(q, f)

        # @show "With sums"
        # @btime lbm.momentum($q, $f)
        # @show "Without allocations"
        j = [0.0, 0.0]
        # @time lbm.momentum!(q, f, j)
        # @btime lbm.momentum!($q, $f, $j)
        # @btime lbm.momentum2!($q, $f, $j)
        lbm.momentum!(q, f, j)
        @test all(j .≈ u)

        # @show lbm.temperature(q, f, ρ, u) * q.speed_of_sound_squared / 2.0
        # @test all(lbm.temperature(q, f, ρ, u) ./ ρ .≈ T)
        # @show lbm.total_energy(q, f)
        # @show lbm.internal_energy(q, f, ρ, u)
        # @show lbm.kinetic_energy(q, f, ρ, u)


    end

    @testset "SRT collision operator" begin
        # q = D2Q9()
        N = 1
        ρ = fill(1.0, N, N)
        u = fill(0.01, N, N, lbm.dimension(q))
        T = fill(0.3, N, N) #./ q.speed_of_sound_squared
        # When τ = 1.0 we immediatly get the equilibrium distribution
        τ = 1.0
        collision_model = SRT(τ)

        f = lbm.equilibrium(q, ρ, u, T);
        # @show q.speed_of_sound_squared
        # @show 1 / q.speed_of_sound_squared
        # @show ρ_1 = lbm.density(q, f)
        # @show j_1 = lbm.momentum(q, f)[:]
        # @show lbm.temperature(q, f, ρ_1, j_1 ./ ρ_1)

        # @show "COLLIDE"
        f_col = collide(collision_model, q, f)

        # @show ρ = lbm.density(q, f)
        # @show j = lbm.momentum(q, f)
        # @show lbm.temperature(q, f, ρ, j ./ ρ)
        @show q
        @test all(f .≈ f_col)
        @inferred lbm.equilibrium(q, ρ, u, T);
        @inferred lbm.equilibrium(q, ρ, u, 1.0);
        # @code_warntype lbm.equilibrium(q, ρ, u, 1.0);
        # @code_warntype collide(collision_model, q, f)
        @inferred collide(collision_model, q, f)
    end

    @testset "Stream and collide" begin
        N = 3
        ρ = fill(1.0, N, N)
        u = fill(0.1, N, N, lbm.dimension(q))
        T = fill(1.0, N, N)
        # When τ = 1.0 we immediatly get the equilibrium distribution
        τ = 1.0
        collision_model = SRT(τ)

        f = lbm.equilibrium(q, ρ, u, T);

        f_in = copy(f)
        for t = 0 : 3
            f_out = collide(SRT(τ), q, f_in)
            f_in = stream(q, f_out)
        end

        @test all(f .≈ f_in)
    end

    @testset "Equilibrium performance" begin
        N = 3
        ρ = fill(1.0, N, N)
        u = fill(0.1, N, N, lbm.dimension(q))
        T = fill(1.0, N, N)
        # When τ = 1.0 we immediatly get the equilibrium distribution

        # @show q
        # @btime lbm.equilibrium($q, $ρ, $u, $T)
        f_1 = lbm.equilibrium(q, ρ, u, T);

        f_2 = Array{Float64, 3}(undef, size(ρ,1), size(ρ,2), length(q.weights));
        # @btime lbm.equilibrium!($q, $ρ, $u, $T, $f_2);
        lbm.equilibrium!(q, ρ, u, T, f_2);

        @test all(f_1 .≈ f_2)
    end

    @testset "Collision performance" begin
        N = 100
        N = 1
        ρ = fill(1.0, N, N)
        u = fill(0.1, N, N, lbm.dimension(q))
        T = fill(1.0, N, N)

        τ = 1.0
        collision_model = SRT(τ)
        # When τ = 1.0 we immediatly get the equilibrium distribution
        f = lbm.equilibrium(q, ρ, u, T);
        f_in = copy(f)

        f_out = copy(f_in)
        f_out_2 = copy(f_in)

        collision_model = SRT(τ)
        # @show q
        # @show "Collide"
        # @btime $f_out = collide($collision_model, $q, $f_in)
        # @show "Collide 2"
        # @btime lbm.collide_2($collision_model, $q, $f_in, $f_out_2)
        # @show to
        @time lbm.collide_2(SRT(τ), q, f_in, f_out_2)

        feq = Array{Float64, 3}(undef, size(f_in,1), size(f_in,2), length(q.weights));
        f_out_3 = copy(f_in)
        # @show "Collide 3"
        # @btime lbm.collide_3($collision_model, $q, $f_in, $f_out_3, $feq)

        f_out = collide(collision_model, q, f_in)
        lbm.collide_2(collision_model, q, f_in, f_out_2)
        @test all(f_out .≈ f_out_2)
    end
end

@testset "Opposite directions of $q" for q in quadratures
    for idx = 1:length(q.weights)
        opposite_idx = opposite(q, idx)

        @test all(q.abscissae[:, idx] .+ q.abscissae[:, opposite_idx] .== 0)
    end
end
