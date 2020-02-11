@testset "Poiseuille Problem" begin #for q in quadratures
    q = D2Q9()

    scale = 1
    ν = 1.0 / 6.0
    problem = PoiseuilleFlow(ν, scale, static = true)
    f_in = LatticeBoltzmann.initialize(AnalyticalEquilibrium(), q, problem)
    f_out = copy(f_in)
    collision_operator = CollisionModel(SRT, q, problem)

    feq = Array{Float64}(undef, size(f_in, 3))
    f = Array{Float64}(undef, size(f_in, 3))
    u = zeros(dimension(q))


    for x_idx in 2, y_idx = 1:size(f_in, 2)
        f = f_in[x_idx, y_idx, :]

        ρ = LatticeBoltzmann.density(q, f)
        LatticeBoltzmann.velocity!(q, f, ρ, u)
        T = LatticeBoltzmann.temperature(q, f, ρ, u)
    end

    t = 0
    Δt = LatticeBoltzmann.delta_t(problem)
    o = size(f_in, 2)
    p = size(f_in, 2) - 1
    LatticeBoltzmann.collide!(collision_operator, q, time = t * Δt, f_old = f_in, f_new = f_out)
    LatticeBoltzmann.stream!(q, f_new = f_in, f_old = f_out)

    for x_idx in 2, y_idx = 1:size(f_in, 2)
        f = f_in[x_idx, y_idx, :]

        ρ = LatticeBoltzmann.density(q, f)
        LatticeBoltzmann.velocity!(q, f, ρ, u)
        T = temperature(q, f, ρ, u)
    end

    for x_idx in 2, y_idx = 1:size(f_in, 2)
        f = f_in[x_idx, y_idx, :]

        ρ = LatticeBoltzmann.density(q, f)
        LatticeBoltzmann.velocity!(q, f, ρ, u)
        T = temperature(q, f, ρ, u)
    end

    # TODO: when an a thermal model is considered, the density fluctuates
    # TODO: when a thermal model is considerd we have constant density
end
