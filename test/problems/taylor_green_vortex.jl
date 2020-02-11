import LinearAlgebra: I, tr

@testset "Taylor Green Vortex Problem" begin #for q in quadratures
    q = D2Q9()

    scale = 1
    ν = 1.0 / 6.0
    problem = TaylorGreenVortex(ν, scale)
    f = LatticeBoltzmann.initialize(AnalyticalEquilibrium(), q, problem)

    # TODO: when an a thermal model is considered, the density fluctuates
    # TODO: when a thermal model is considerd we have constant density
    @testset "Hermite Pressure component" begin
        q = D2Q9()
        D = dimension(q)
        N = div(LatticeBoltzmann.order(q), 2)
        Hs = [
            [hermite(Val{n}, q.abscissae[:, i], q) for i = 1:length(q.weights)] for n = 1:N
        ]

        scale = 1
        ν = 1.0 / 6.0
        problem = DecayingShearFlow(ν, scale, static = false)
        problem = TaylorGreenVortex(ν, scale, static = true)

        x = 0.30
        y = 0.30
        f = LatticeBoltzmann.initial_condition(AnalyticalEquilibrium(), q, problem, x, y)

        ρ = sum(f)

        τ = problem.ν
        τ = q.speed_of_sound_squared * LatticeBoltzmann.lattice_viscosity(problem) + 0.5

        a_f = [sum([f[idx] * Hs[n][idx] for idx = 1:length(q.weights)]) for n = 1:N]

        P = a_f[2] - (a_f[1] * a_f[1]') / ρ - ρ * I
        P = P * (1 - 1 / (2 * τ))

        σ_lb = (P - I * tr(P) / (D)) / (problem.u_max)
        σ = LatticeBoltzmann.deviatoric_tensor(q, problem, x, y, 0.0)
    end
    @testset "Pressure component" begin
        scale = 1
        ν = 2.0 / 6.0
        problem = TaylorGreenVortex(ν, scale, static = false)
        problem = DecayingShearFlow(ν, scale, static = false)

        x = 0.30
        y = 0.30
        f = LatticeBoltzmann.initial_condition(AnalyticalEquilibrium(), q, problem, x, y)

        ρ = LatticeBoltzmann.density(q, f)

        cs = 1 / q.speed_of_sound_squared
        # cs = q.speed_of_sound_squared
        d_u = problem.u_max * LatticeBoltzmann.velocity_gradient(problem, x, y, 0.0)
        τ = problem.ν
        τ = q.speed_of_sound_squared * LatticeBoltzmann.lattice_viscosity(problem) + 0.5
        for f_idx = 1:length(f)
            f[f_idx] +=
                -(q.weights[f_idx] * ρ * τ / cs) *
                sum(LatticeBoltzmann.hermite(Val{2}, q.abscissae[:, f_idx], q) .* d_u)
        end


        u = zeros(LatticeBoltzmann.dimension(q))
        LatticeBoltzmann.velocity!(q, f, ρ, u)
        T = LatticeBoltzmann.temperature(q, f, ρ, u)

        σ_lb = LatticeBoltzmann.momentum_flux(q, f, ρ, u) - I * LatticeBoltzmann.pressure(q, f, ρ, u)
        σ = LatticeBoltzmann.deviatoric_tensor(q, problem, x, y, 0.0)
    end
end
