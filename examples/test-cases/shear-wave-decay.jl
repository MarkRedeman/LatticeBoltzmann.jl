using lbm, Plots, DataFrames
# using PGFPlotsX

import lbm: StopCriteria, CompareWithAnalyticalSolution, TrackHydrodynamicErrors

# ProcessingMethod(problem::DecayingShearFlow, should_process, n_steps, stop_criteria = StopCriteria(problem)) =
#     CompareWithAnalyticalSolution(problem, should_process, n_steps, stop_criteria)
# ProcessingMethod(problem::DecayingShearFlow, should_process, n_steps, stop_criteria = StopCriteria(problem)) =
#     TrackHydrodynamicErrors(problem, should_process, n_steps, stop_criteria)

function shear_wave_convergence_analysis()
    q = D2Q9()
    ν = 2.0 / 6.0

    stats = DataFrame(
        [Float64[], Int[], Any[]],
        [:ν, :scale, :stats]
    )

    # for scale = [1//2, 1, 2, 4]
    #     problem = DecayingShearFlow(ν, scale, static = true)

    #     Δt = delta_t(problem)
    #     t_end = 1.0
    #     n_steps = round(Int, t_end / Δt)

    #     # CompareWithAnalyticalSolution
    #     process_method = TrackHydrodynamicErrors(
    #         problem,
    #         false,
    #         n_steps,
    #         StopCriteria(problem)
    #     )


    #     result = lbm.simulate(problem, q, process_method = process_method, t_end = t_end)

    #     push!(stats, [ν, scale, result.processing_method.df[end]])
    # end

    # scales = 1:4
    scales = [1//2, 1, 2, 4]
    results = map(scales) do scale
        problem = DecayingShearFlow(ν, scale, static = false)
        τ = 1.0
        problem = lbm.TGV(q, τ, scale)
        u_max = 0.01 / scale
        # ν = τ / (2.0 * q.speed_of_sound_squared)
        NX = Int(scale * 8)
        NY = Int(scale * 8)

        problem = lbm.TGV(q, τ, scale, NX, NY, u_max)


        Δt = delta_t(problem)
        t_end = 1.0
        t_end = 0.250
        t_end = 100
        n_steps = round(Int, t_end / Δt)

        # CompareWithAnalyticalSolution
        process_method = TrackHydrodynamicErrors(
            problem,
            false,
            n_steps,

            StopCriteria(problem)
        )

        res = lbm.simulate(problem, q, process_method = process_method, t_end = t_end)
        @show res.processing_method.df[end]

        return res.processing_method.df[end]
    end

    @show results

    p = (
        plot(map(scale -> 8*scale, scales), getfield.(results, :error_u), scale=:log10, title="u"),
        plot(map(scale -> 8*scale, scales), getfield.(results, :error_σ_xx), scale=:log10, title="sigma xx"),
        plot(map(scale -> 8*scale, scales), getfield.(results, :error_σ_xy), scale=:log10, title="sigma xy"),
        plot(map(scale -> 8*scale, scales), getfield.(results, :error_p), scale=:log10, title="p"),
    )

    p = plot(map(scale -> 8*scale, scales), getfield.(results, :error_u), scale=:log10, label="u")
    # plot!(p, map(scale -> 8*scale, scales), getfield.(results, :error_σ_xx), scale=:log10, title="sigma xx")
    plot!(p, map(scale -> 8*scale, scales), getfield.(results, :error_σ_xy), scale=:log10, label="sigma xy")
    plot!(p, map(scale -> 8*scale, scales), getfield.(results, :error_p), scale=:log10, label="p")
    plot!(p, x -> 0.1 * x.^(-2), scale=:log10, label="x^-2")
    plot!(p, x -> 0.1 * x.^(-4), scale=:log10, label="x^-4")

    return results, p
end

function shear_wave_optimal_relaxation_time()
end

function shear_wave_optimal_two_relaxation_time()
end

function shear_wave_multispeed_quadrature_analysis()
end

result = shear_wave_convergence_analysis()
