module ExamplesShearWaveDecay

include("lbm.jl")
# using PGFPlotsX

# using PGFPlots
# using PGFPlotsX
# pgfplots() # x?

# ProcessingMethod(problem::DecayingShearFlow, should_process, n_steps, stop_criteria = StopCriteria(problem)) =
#     CompareWithAnalyticalSolution(problem, should_process, n_steps, stop_criteria)
# ProcessingMethod(problem::DecayingShearFlow, should_process, n_steps, stop_criteria = StopCriteria(problem)) =
#     TrackHydrodynamicErrors(problem, should_process, n_steps, stop_criteria)

# function shear_wave_convergence_analysis()
#     q = D2Q9()
#     ν = 2.0 / 6.0

#     stats = DataFrame(
#         [Float64[], Int[], Any[]],
#         [:ν, :scale, :stats]
#     )

#     # for scale = [1//2, 1, 2, 4]
#     #     problem = DecayingShearFlow(ν, scale, static = false)

#     #     Δt = delta_t(problem)
#     #     t_end = 1.0
#     #     n_steps = round(Int, t_end / Δt)

#     #     # CompareWithAnalyticalSolution
#     #     process_method = TrackHydrodynamicErrors(
#     #         problem,
#     #         false,
#     #         n_steps,
#     #         StopCriteria(problem)
#     #     )


#     #     result = lbm.simulate(problem, q, process_method = process_method, t_end = t_end)

#     #     push!(stats, [ν, scale, result.processing_method.df[end]])
#     # end

#     # scales = 1:4
#     scales = [1//2, 1, 2, 4]
#     results = map(scales) do scale
#         problem = DecayingShearFlow(ν, scale, static = false)
#         τ = 1.0
#         problem = lbm.TGV(q, τ, scale)
#         u_max = 0.01 / scale
#         # ν = τ / (2.0 * q.speed_of_sound_squared)
#         NX = Int(scale * 8)
#         NY = Int(scale * 8)

#         problem = lbm.TGV(q, τ, scale, NX, NY, u_max)


#         Δt = delta_t(problem)
#         t_end = 1.0
#         t_end = 0.250
#         t_end = 100
#         n_steps = round(Int, t_end / Δt)

#         # CompareWithAnalyticalSolution
#         process_method = TrackHydrodynamicErrors(
#             problem,
#             false,
#             n_steps,

#             StopCriteria(problem)
#         )

#         res = lbm.simulate(problem, q, process_method = process_method, t_end = t_end)
#         @show res.processing_method.df[end]

#         return res.processing_method.df[end]
#     end

#     @show results

#     p = (
#         plot(map(scale -> 8*scale, scales), getfield.(results, :error_u), scale=:log10, title="u"),
#         plot(map(scale -> 8*scale, scales), getfield.(results, :error_σ_xx), scale=:log10, title="sigma xx"),
#         plot(map(scale -> 8*scale, scales), getfield.(results, :error_σ_xy), scale=:log10, title="sigma xy"),
#         plot(map(scale -> 8*scale, scales), getfield.(results, :error_p), scale=:log10, title="p"),
#     )

#     p = plot(map(scale -> 8*scale, scales), getfield.(results, :error_u), scale=:log10, label="u")
#     # plot!(p, map(scale -> 8*scale, scales), getfield.(results, :error_σ_xx), scale=:log10, title="sigma xx")
#     plot!(p, map(scale -> 8*scale, scales), getfield.(results, :error_σ_xy), scale=:log10, label="sigma xy")
#     plot!(p, map(scale -> 8*scale, scales), getfield.(results, :error_p), scale=:log10, label="p")
#     plot!(p, x -> 0.1 * x.^(-2), scale=:log10, label="x^-2")
#     plot!(p, x -> 0.1 * x.^(-4), scale=:log10, label="x^-4")

#     return results, p
# end
# logspace(a, b, length) = 10.^(range(log10(a), stop = log10(b), length = length))

function shear_wave_optimal_relaxation_time()
    # τ = range()0.5:1.0

    let
        s3 = []
        # for q in (D2Q9 = D2Q9(), D2Q13 = D2Q13(), D2Q17 = D2Q17(), D2Q21 = D2Q21(), D2Q37 = D2Q37())
        for q in [D2Q9()]
            # for q in lbm.Quadratures
            stats = DataFrame([Float64[], Float64[], Float64[], Quadrature[]], [:τ, :ν, :error_u, :q])
            # for τ = 0.5:0.01:100.0
            τs = 0.5 .+ [10^i for i in range(-8, 1.5, length = 1000)]
            for τ = τs
                ν = (τ - 0.5) / q.speed_of_sound_squared

                # problem = DecayingShearFlow(ν, 2, static = true)
                # problem = PoiseuilleFlow(ν, 2, static = true)
                problem = CouetteFlow(ν, 2, static = true)

                result = lbm.simulate(
                    problem,
                    q,
                    t_end = 1.0,
                    should_process = false,
                    collision_model=SRT
                )
                if (! isnan(result.processing_method.df[end].error_u))
                    push!(stats, [τ, ν, result.processing_method.df[end].error_u, q])
                end
            end

            push!(s3, stats)
        end

        p2 = plot(xlabel = L"\tau", ylabel = L"\epsilon_{u}", scale=:log10, legend=:bottomright, title = L"\gau");
        for stats in s3
            min_idx = argmin(stats.error_u)
            @show stats.τ[min_idx]
            @show stats.error_u[min_idx]
            plot!(p2, stats.τ, stats.error_u, label = string(stats.q[1]))
            @show stats.τ[min_idx], stats.error_u[min_idx]
            scatter!(p2, [stats.τ[min_idx]], [stats.error_u[min_idx]], label = string(stats.q[1]), markercolor = :white, markeralpha = 0.3)
        end
        display(p2)

        p1 = plot(xlabel = L"\nu", ylabel = L"\epsilon_{u}", legend=:bottomright, ylim=(-Inf, 100), title = L"\nu");
        for stats in s3
            min_idx = argmin(stats.error_u)
            @show stats.τ[min_idx]
            @show stats.error_u[min_idx]
            plot!(p1, stats.ν, stats.error_u, label = string(stats.q[1]))
            scatter!(p1, [stats.ν[min_idx]], [stats.error_u[min_idx]], label = string(stats.q[1]))
        end
        display(p1)
        return p, s3
    end
end

function shear_wave_optimal_two_relaxation_time()
end

function shear_wave_multispeed_quadrature_analysis()
end

function shear_wave_velocity_profile(
    q = D2Q9(),
    initialization_strategy = nothing,
    τ = 1.0,
    scale = 2
)
    ν = τ / (2.0 * q.speed_of_sound_squared)
    scale = 1
    problem = DecayingShearFlow(ν, scale, static = false)

    if isnothing(initialization_strategy)
        initialization_strategy = InitializationStrategy(problem)
    end

    @show q initialization_strategy

    t_end = 1.20
    t_end = 1.0
    t_end = 2.250
    plot_every = 10000

    simulate(
        problem,
        q,
        # should_process = false,
        initialization_strategy = initialization_strategy,
        t_end = 1.0,
        # process_method = ShowVelocityError(problem, plot_every, Float64[]),
    )
end

function shear_wave_convergence_analysis(q = D2Q9(), initialization_strategy = AnalyticalEquilibrium(), τ = 1.0; scales = [1, 2, 4, 8])
    ν = τ / (2.0 * q.speed_of_sound_squared)

    t_end = 0.20
    t_end = 1.0

    results = map(scales) do scale
        problem = DecayingShearFlow(ν, scale, static = true)

        Δt = delta_t(problem)
        n_steps = round(Int, t_end / Δt)

        # CompareWithAnalyticalSolution
        process_method = TrackHydrodynamicErrors(
            problem,
            false,
            n_steps,
            lbm.NoStoppingCriteria()
        )
        @show scale problem initialization_strategy n_steps

        @time res = simulate(
            problem,
            q,
            process_method = process_method,
            initialization_strategy = initialization_strategy,
            t_end = t_end,
        )

        return res.processing_method.df[end]
    end

    return (
        quadrature = q,
        scales = scales,
        results = results
    )
end

function main(τ = 1.0, scale = 2)
    quadratures = [
        # D2Q4(),
        # D2Q5(),
        D2Q9(),
        D2Q13(),
        D2Q17(),
        D2Q21(),
        D2Q37(),
    ]

    iteration_strategies = [
        "Constant density" => ConstantDensity(),
        "Iterative" => IterativeInitializationMeiEtAl(1.0τ, 1E-11),
        "Analytical equilibrum" => AnalyticalEquilibrium(),
        "Analytical offequilibrium" => AnalyticalEquilibriumAndOffEquilibrium(),

        # # NOTE: requires a static (forced) TGV
        # "Zero initial velocity" => ZeroVelocityInitialCondition(),
    ]
    plots = []
    for itt in iteration_strategies
        results = map(quadratures) do q
            shear_wave_velocity_profile(
                q,
                # AnalyticalEquilibrium(),
                itt[2],
                τ,
                scale
            )
        end
        p = plot_error_progresion(results)
        plot!(p, title = itt[1])
        push!(plots, p)
    end
    display(plot(plots..., legend = :bottomleft, size=(900, 600)))
    results = map(quadratures) do q
        shear_wave_velocity_profile(
            q,
            AnalyticalEquilibrium(),
            τ,
            scale
        )
    end

    # scales = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    # scales = [1, 2, 4, 8, 12, 16]
    scales = [1, 2, 4, 8, 16, 32, 64]

    scales = [1, 2, 4, 8, 16]
    scales = [1, 2, 4, 8, 16, 32]
    # scales = [1//2,1, 2, 4, 8, 16]

    scales = [1, 2, 4, 8]
    iteration_strategies = [
        ConstantDensity(),
        AnalyticalEquilibrium(),
        AnalyticalEquilibriumAndOffEquilibrium(),
        IterativeInitializationMeiEtAl(τ, 1E-11),

        # NOTE: requires a static (forced) TGV
        ZeroVelocityInitialCondition(),
    ]
    # iteration_strategies = [
    #     "Constant density" => ConstantDensity(),
    #     "Iterative" => IterativeInitializationMeiEtAl(1.0τ, 1E-7),
    #     "Analytical equilibrum" => AnalyticalEquilibrium(),
    #     "Analytical offequilibrium" => AnalyticalEquilibriumAndOffEquilibrium(),

    #     # # NOTE: requires a static (forced) TGV
    #     # "Zero initial velocity" => ZeroVelocityInitialCondition(),
    # ]

    convergence_results = map(quadratures) do q
        shear_wave_convergence_analysis(
            q,
            iteration_strategies[1],
            τ,
            scales = scales
        )
    end
    convergence_results_iterative = map(quadratures) do q
        shear_wave_convergence_analysis(
            q,
            iteration_strategies[4],
            τ,
            scales = scales
        )
    end
    convergence_results_equilibrium = map(quadratures) do q
        shear_wave_convergence_analysis(
            q,
            iteration_strategies[2],
            τ,
            scales = scales
        )
    end
    convergence_results_offequilibrium = map(quadratures) do q
        shear_wave_convergence_analysis(
            q,
            iteration_strategies[3],
            τ,
            scales = scales
        )
    end

    return (
        results = results,
        convergence_results = convergence_results,
        convergence_results_iterative = convergence_results_iterative,
        convergence_results_equilibrium = convergence_results_equilibrium,
        convergence_results_offequilibrium = convergence_results_offequilibrium,
    )
end

function plot_main(main = main())
    error_location = plot_error_locations(main.results)
    error_progression = plot_error_progresion(main.results)

    # display(error_location)
    display(error_progression)

    for plot_moment = [:error_u, :error_p, :error_σ_xx, :error_σ_xy]
        convergence = plot_convergence(main.convergence_results, plot_moment, 8)
        convergence_iterative = plot_convergence(main.convergence_results_iterative, plot_moment, 8)
        convergence_equilibrium = plot_convergence(main.convergence_results_equilibrium, plot_moment, 8)
        convergence_offequilibrium = plot_convergence(main.convergence_results_offequilibrium, plot_moment, 8)
        display(
            plot(
                plot!(convergence, title = "Constant density"),
                plot!(convergence_equilibrium, title = "Analytical equilibrium"),
                plot!(convergence_offequilibrium, title = "Analytical offequilibrium"),
                plot!(convergence_iterative, title = "Iterative"),
                legend = nothing,
                size = (900, 600)
            )
        )
    end

    convergence = plot_convergence(main.convergence_results, :error_u, 8)
    convergence_iterative = plot_convergence(main.convergence_results_iterative, :error_u, 8)
    convergence_equilibrium = plot_convergence(main.convergence_results_equilibrium, :error_u, 8)
    convergence_offequilibrium = plot_convergence(main.convergence_results_offequilibrium, :error_u, 8)


    return (
        results = main.results,
        convergence_results = main.convergence_results,
        convergence_results_iterative = main.convergence_results_iterative,
        convergence_results_equilibrium = main.convergence_results_equilibrium,
        convergence_results_offequilibrium = main.convergence_results_offequilibrium,

        plots = (
            error_location = error_location,
            error_progression = error_progression,
            convergence = convergence,
        )
    )
end

# result = shear_wave_convergence_analysis()
#

function thesis_main(q = D2Q9(), scale = 1, initialization_strategy = AnalyticalEquilibrium())
    # 1.0 Show solutions at time steps (for a given quadrature, viscosity)
    # t = [0.0, 0.25, 0.5, 0.75, 1.0]
    static = false

    snapshots_at = [0.0, 0.25, 0.5, 0.75, 1.0]
    snapshots_at = [0.0, 0.05, 0.15, 0.25]
    ν = 1.0 / (2 * q.speed_of_sound_squared)

    problem = DecayingShearFlow(ν, scale, static = false, A = 3.0)
    Δt = delta_t(problem)

    snapshot_results = map(snapshots_at) do t_end
        simulate(problem, q, should_process = false, initialization_strategy = initialization_strategy, t_end = t_end,)
    end

    snapshot_plots = plot_snapshots(problem, snapshot_results, snapshots_at, q)

    # 2.0 For the same configuration show convergence rate of
    # velocity, pressure, stress tensor
    convergence_plots = nothing

    # 3.0a Repeat 1 - 2 for a static version (to check the results)
    # 3.0b Find optimal relaxation time for SRT

    # 4.0 Find optimal relaxation times for TRT

    # 5.0 Possibly check for other quadratures?

    # 6.0 Check effect of u_max vs τ / ν

    return (
        snapshot_plots = snapshot_plots,
        convergence_plots = convergence_plots
    )
end

function plot_snapshots(problem, snapshot_results, snapshots, q)
    velocity_profile_x = plot(xlabel = "x", ylabel=L"u_x")
    velocity_profile_y = plot(xlabel = "x", ylabel=L"u_y", legend=:bottomright, title = string("Velocity profile at ", latexstring("y = \\pi")))
    pressure_profile = plot(xlabel = "x", ylabel=L"p")
    temperature_profile = plot(xlabel = "x", ylabel=L"T")
    sigma_xx_profile = plot(xlabel = "x", ylabel=L"\sigma_{xx}")
    sigma_xy_profile = plot(xlabel = "x", ylabel=L"\sigma_{xy}")

    for (time, snapshot) in zip(snapshots, snapshot_results)
        f_in = snapshot.f_stream

        # Pre-allocate Macroscopic Variables
        ρ = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))
        u = Array{Float64}(undef, size(f_in, 1), size(f_in, 2), dimension(q))
        p = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))
        T = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))
        σ_xx = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))
        σ_xy = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))

        Nx = size(f_in, 1)
        Ny = size(f_in, 2)

        x_range, y_range = range(problem)

        f = Array{Float64}(undef, size(f_in, 3))
        u_ = zeros(dimension(q))
        @inbounds for x_idx = 1:Nx, y_idx = 1:Ny
            @inbounds for f_idx = 1:size(f_in, 3)
                f[f_idx] = f_in[x_idx, y_idx, f_idx]
            end
            x = x_range[x_idx]
            y = y_range[y_idx]

            ρ_ = density(q, f)
            velocity!(q, f, ρ_, u_)
            T_ = temperature(q, f, ρ_, u_)
            p_ = pressure(q, f, ρ_, u_)

            τ = q.speed_of_sound_squared * lbm.lattice_viscosity(problem)
            σ_ = lbm.deviatoric_tensor(q, τ, f, ρ_, u_)

            ρ_ = lbm.dimensionless_density(problem, ρ_)
            u_ = lbm.dimensionless_velocity(problem, u_)
            T_ = lbm.dimensionless_temperature(q, problem, T_)
            p_ = lbm.dimensionless_pressure(q, problem, p_)
            ρ[x_idx, y_idx] = ρ_
            u[x_idx, y_idx, :] = u_
            p[x_idx, y_idx] = p_
            T[x_idx, y_idx] = T_

            σ_ = lbm.dimensionless_stress(problem, σ_)
            σ_xx[x_idx, y_idx] = σ_[1, 1]
            σ_xy[x_idx, y_idx] = σ_[1, 2]
        end
        s = (1000, 500)
        y_pos = round(Int, problem.NY / 2)
        domain = x_range[1:Nx]


        x_range, y_range = range(problem)
        velocity = (x, y, t) -> lbm.velocity(problem, x, y, t)
        σ = (x, y, t) -> lbm.deviatoric_tensor(q, problem, x, y, t)
        pr = (x, y, t) -> lbm.pressure(q, problem, x, y, t)

        exact_range = range(0.0, length = 1000, stop = problem.domain_size[1])

        scatter!(velocity_profile_x, domain, u[:, y_pos, 1], label = latexstring("t = ", time))
        plot!((x) -> velocity(x, y_range[y_pos], time)[1], exact_range, label = "", linecolor = :gray, linealpha = 0.2, linestyle = :dash)

        scatter!(velocity_profile_y, domain, u[:, y_pos, 2], label = latexstring("t = ", time), markershape = :auto, markersize = 6)
        plot!((x) -> velocity(x, y_range[y_pos], time)[2], exact_range, label = "", linecolor = :gray, linealpha = 0.2, linestyle = :dash)

        scatter!(pressure_profile, domain, p[:, y_pos], label = latexstring("t = ", time))
        plot!(pressure_profile, (x) -> pr(x, y_range[y_pos], time), exact_range, label = "", linecolor = :gray, linealpha = 0.2, linestyle = :dash)

        plot!(temperature_profile, domain, T[:, y_pos], label = latexstring("t = ", time))

        scatter!(sigma_xx_profile, domain, σ_xx[:, y_pos], label = latexstring("t = ", time))
        plot!(sigma_xx_profile, (x) -> σ(x, y_range[y_pos], time)[1,1], exact_range, label = "", linecolor = :gray, linealpha = 0.2, linestyle = :dash)

        scatter!(sigma_xy_profile, domain, σ_xy[:, y_pos], label = latexstring("t = ", time))
        plot!(sigma_xy_profile, (x) -> σ(x, y_range[y_pos], time)[1,2], exact_range, label = "", linecolor = :gray, linealpha = 0.2, linestyle = :dash)
    end

    return (
        velocity_profile_x = velocity_profile_x,
        velocity_profile_y = velocity_profile_y,
        temperature_profile = temperature_profile,
        pressure_profile = pressure_profile,
        sigma_xx_profile = sigma_xx_profile,
        sigma_xy_profile = sigma_xy_profile,
    )
end

function show_velocity_evolution_of_steady_state()
    problem = DecayingShearFlow(1.0 / 6.0, 16, A = 0.5, static = true)
    ν = lbm.viscosity(problem)
    Δt = delta_t(problem);
    snapshot_at = round.(Int, [
        0.01 * 1.0 / (ν * Δt),
        # 0.05 * 1.0 / (ν * Δt),
        0.1 * 1.0 / (ν * Δt),
        # 0.5 * 1.0 / (ν * Δt),
        1.0 * 1.0 / (ν * Δt),
        # 5.0 / (ν * Δt),
        10.0 / (ν * Δt)
    ])

    q = D2Q9()
    result = simulate(
        problem,
        q,
        t_end = 10.0 / lbm.viscosity(problem),
        process_method = lbm.TakeSnapshots(problem, snapshot_at),
        initialization_strategy = lbm.ZeroVelocityInitialCondition()
    );

    p = lbm.visualize(result.processing_method, q)
    plot!(p.velocity_profile_x, legend = nothing)
    return p
end

end
