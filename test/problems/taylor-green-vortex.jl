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

        @show a_f[2]
        P = a_f[2] - (a_f[1] * a_f[1]') / ρ - ρ * I
        P = P * (1 - 1 / (2 * τ))

        σ_lb = (P - I * tr(P) / (D)) / (problem.u_max)
        σ = LatticeBoltzmann.deviatoric_tensor(q, problem, x, y, 0.0)

        @show σ σ_lb
        @show σ - σ_lb
        @show σ ./ σ_lb
        @show σ_lb ./ σ
        # return


        # @inbounds for f_idx = 1 : nf
        #     f_out[x_idx, y_idx, f_idx] = q.weights[f_idx] * (
        #         ρ +
        #         cs * ρ * sum(u .* Hs[1][f_idx]) +
        #         sum([
        #             cs^n * dot(a_coll[n], Hs[n][f_idx]) / (factorial(n))
        #             for n = 2:N
        #         ])
        #     )
        # end
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
        @show d_u
        τ = q.speed_of_sound_squared * LatticeBoltzmann.lattice_viscosity(problem) + 0.5
        @show τ
        for f_idx = 1:length(f)
            f[f_idx] +=
                -(q.weights[f_idx] * ρ * τ / cs) *
                sum(LatticeBoltzmann.hermite(Val{2}, q.abscissae[:, f_idx], q) .* d_u)
        end
        @show 1 / delta_t(problem)
        @show 1 / LatticeBoltzmann.delta_x(problem)
        @show 1 / viscosity(problem)


        u = zeros(LatticeBoltzmann.dimension(q))
        LatticeBoltzmann.velocity!(q, f, ρ, u)
        T = LatticeBoltzmann.temperature(q, f, ρ, u)
        @show LatticeBoltzmann.momentum_flux(q, f, ρ, u)
        @show LatticeBoltzmann.pressure(q, f, ρ, u)

        @show -(q.speed_of_sound_squared / problem.ν) *
              (LatticeBoltzmann.momentum_flux(q, f, ρ, u) - I * LatticeBoltzmann.pressure(q, f, ρ, u)) ./
              problem.u_max
        @show ρ * u * u'

        @show LatticeBoltzmann.deviatoric_tensor(q, problem, x, y, 0.0)
        # @show LatticeBoltzmann.pressure_tensor(q, problem, x, y, 0.0)
        @show problem.u_max
        @show LatticeBoltzmann.pressure_tensor(q, problem, x, y, 0.0) -
              LatticeBoltzmann.deviatoric_tensor(q, problem, x, y, 0.0)
        @show LatticeBoltzmann.pressure(q, f, ρ, u) - LatticeBoltzmann.pressure(q, problem, x, y, 0.0)


        @warn "Check deviatoric stress"
        σ_lb = LatticeBoltzmann.momentum_flux(q, f, ρ, u) - I * LatticeBoltzmann.pressure(q, f, ρ, u)
        σ = LatticeBoltzmann.deviatoric_tensor(q, problem, x, y, 0.0)
        @show σ_lb σ
        @show (1 - 1 / (2 * τ))
        @show 0.5 * (1 / cs) * (1 / problem.u_max) * σ_lb / (1 - 1 / (2 * τ))
        # @show (1 / (problem.ν^2 * problem.u_max)) * σ_lb / (1 + 1 / (2 * problem.ν))
        # @show (q.speed_of_sound_squared / (problem.ν * problem.u_max) ) * (
        #     σ_lb
        # )
        # @show σ

        # @show (1 / (problem.ν^2 * problem.u_max)) * σ_lb / (1 + 1 / (2 * problem.ν)) + σ/2
        # @show (q.speed_of_sound_squared / (problem.ν * problem.u_max) ) * (
        #     σ_lb
        # ) + σ
        # return

    end
end
