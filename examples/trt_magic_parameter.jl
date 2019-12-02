using DataFrames
using lbm
using Plots

# let
    q = D2Q9()
    # problem = PoiseuilleFlow(ν, 2, static = true)
    s = []

    scale_range = [1, 2, 4]
    scale_range = [1, 2]
    scale_range = [4]
    scale_range = [2]
    Λ_range = vcat(
        range(0.01, stop = 0.75, step = 0.0005),
        range(0.75, stop = 1.0, step = 0.05)
    )

    Quadratures = lbm.Quadratures
    Quadratures = [D2Q9()]

    for q in Quadratures
    for τ = 0.5:0.01:2.0
    ν = τ / (2 * q.speed_of_sound_squared)
    # ν = 1.00 / (2 * q.speed_of_sound_squared)
    for scale = scale_range
        # problem = DecayingShearFlow(ν, scale, static = false)
        problem = PoiseuilleFlow(ν, scale, static = true)
        stats = DataFrame([
            Float64[], Float64[], Float64[], Any[], Float64[], Float64[], Float64[]
        ], [
            :Λ, :scale, :τ, :q, :error_u, :error_σ_xx, :error_σ_xy
        ])

        for Λ = Λ_range
            result = lbm.siumlate(
                problem,
                q,
                t_end = 0.15,
                should_process = false,
                collision_model=lbm.TRT_Λ(Λ)
            )

            if (! isnan(result.processing_method.df[end].error_u))
                push!(
                    stats,
                    [
                        Λ,
                        scale,
                        τ,
                        q,
                        result.processing_method.df[end].error_u,
                        result.processing_method.df[end].error_σ_xx,
                        result.processing_method.df[end].error_σ_xy
                    ]
                )
            end
        end

        # plot(stats.Λ, stats.error_u, yscale=:log10)
        # gui()
        push!(s, stats)
    end
    end
    end

    p = plot()
    for stats in s
        @show stats.Λ[argmin(stats.error_u)]
        plot!(stats.Λ, stats.error_u, yscale=:log10, linestyle=:dot)
    end
    gui()


    optimal_Λs = []
    for stats in s
        @show stats.Λ[argmin(stats.error_σ_xy)]
        push!(optimal_Λs, [
              stats.τ[argmin(stats.error_u)],
              stats.Λ[argmin(stats.error_u)],
        ])
    end

    τ_ = 0.5 .+ last.(optimal_Λs) ./ (first.(optimal_Λs) .- 0.5)
    plot(
        first.(optimal_Λs),
        last.(optimal_Λs),
        xlabel = "Tau",
        label="Optimal Lambda"
    )
    plt = twinx()
    plot!(
        plt,
        first.(optimal_Λs),
        0.5 .+ last.(optimal_Λs) ./ (first.(optimal_Λs) .- 0.5),
        xlabel = "Tau",
        label="Optimal Tau -"
    )
    gui()

    # stats.τ[argmin(stats.error_u)] = 1.438
    # 1.438
    # julia> @show stats.τ[argmin(stats.error_u)]
    # stats.τ[argmin(stats.error_u)] = 1.76
    # 1.76
# end
