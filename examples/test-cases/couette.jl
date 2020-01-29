using LatticeBoltzmann, Plots, DataFrames
using LaTeXStrings

module Examples
module Couette

using LatticeBoltzmann, Plots, DataFrames
# using ElectronDisplay
using LaTeXStrings
using JLD2

# using PGFPlots
# using PGFPlotsX
# pgfplots() # x?

import LatticeBoltzmann: StopCriteria,
    CompareWithAnalyticalSolution,
    TrackHydrodynamicErrors,
    ZeroVelocityInitialCondition,
    IterativeInitializationMeiEtAl,
    density,
    velocity!,
    dimensionless_velocity,
    ProcessingMethod,
    next!,
    InitializationStrategy,
    ShowVelocityError

# line_style(q::Quadrature) = (:line, :dot)

line_style(q::D2Q4) = (:line)
line_style(q::D2Q5) = (:line, :dash, 1.0, 3.0)
line_style(q::D2Q9) = (:line)
line_style(q::D2Q13) = (:line)
line_style(q::D2Q17) = (:line)
line_style(q::D2Q21) = (:line)
line_style(q::D2Q37) = (:line, :dash, 1.0, 3.0)
marker_style(q::Quadrature) = ()

marker_style(q::D2Q4) = ()
marker_style(q::D2Q5) = (:hexagon)
marker_style(q::D2Q9) = ()
marker_style(q::D2Q13) = ()
marker_style(q::D2Q17) = ()
marker_style(q::D2Q21) = ()
marker_style(q::D2Q37) = (:hexagon)
# line_style(q::D2Q4) = (linestyle = :dot, linecolor = :red)

function couette_velocity_profile(
    q = D2Q9(),
    initialization_strategy = nothing,
    τ = 1.0,
    scale = 2
)
    ν = τ / (2.0 * q.speed_of_sound_squared)
    problem = CouetteFlow(ν, scale)

    if isnothing(initialization_strategy)
        initialization_strategy = InitializationStrategy(problem)
    end

    @show q initialization_strategy

    t_end = 0.20
    plot_every = 100

    if typeof(initialization_strategy) == ZeroVelocityInitialCondition
        t_end = 10.0
        plot_every = 1000
    end

    simulate(
        problem,
        q,
        should_process = false,
        initialization_strategy = initialization_strategy,
        t_end = t_end,
        process_method = ShowVelocityError(problem, plot_every, Float64[]),
    )
end

function couette_convergence_analysis(q = D2Q9(), initialization_strategy = ZeroVelocityInitialCondition(), τ = 1.0; scales = [1, 2, 4, 8])
    ν = τ / (2.0 * q.speed_of_sound_squared)

    t_end = 0.20
    if typeof(initialization_strategy) == ZeroVelocityInitialCondition
        t_end = 3.0
    end
    t_end = 10.0

    results = map(scales) do scale
        problem = CouetteFlow(ν, scale)

        Δt = delta_t(problem)
        n_steps = round(Int, t_end / Δt)

        # CompareWithAnalyticalSolution
        process_method = TrackHydrodynamicErrors(
            problem,
            false,
            n_steps,
            LatticeBoltzmann.NoStoppingCriteria()
        )
        @show scale problem initialization_strategy n_steps

        @time res = simulate(
            problem,
            q,
            process_method = process_method,
            initialization_strategy = initialization_strategy,
            t_end = t_end
        )

        return res.processing_method.df[end]
    end

    return (
        quadrature = q,
        scales = scales,
        results = results
    )

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


function plot_error_locations(results)
    p = plot()

    for result in results
        f = result.f_stream
        problem = result.processing_method.problem
        q = result.quadrature

        analytical_velocity = (y) -> velocity(problem, 1.0, y)[1]
        x_domain = (0.0, problem.domain_size[1])
        y_domain = (0.0, problem.domain_size[2])

        x_idx = 1
        v_y = zeros(problem.NY)
        v_e = zeros(problem.NY)
        v_a = zeros(problem.NY)
        u = zeros(dimension(q))
        x_range, y_range = range(problem)
        for y_idx = 1 : problem.NY
            ρ = density(q, f[x_idx, y_idx, :])
            velocity!(q, f[x_idx, y_idx, :], ρ, u)
            u = dimensionless_velocity(problem, u)
            v_e[y_idx] = abs(u[1] - analytical_velocity(y_range[y_idx]))
            v_y[y_idx] = u[1]
            v_a[y_idx] = analytical_velocity(y_range[y_idx])
        end

        scatter!(
            p,
            y_range,
            v_e,
            label = string(q),
            marker = marker_style(result.quadrature),
        )
    end

    plot!(
        p,
        legend=:bottomleft,
        ylabel = L"\epsilon_u",
        xlabel = "y",
        yscale = :log10,
    )

    p
end
function plot_error_progresion(results)
    p = plot()

    for result in results
        q = result.quadrature
        problem = result.processing_method.problem
        Δt = delta_t(problem)
        plot!(
            p,
            Δt * (1 : length(result.processing_method.l2_norms)),
            result.processing_method.l2_norms,
            line = line_style(result.quadrature),
            # marker = marker_style(result.quadrature),
            label=LaTeXString(string(q)),
        )
    end

    plot!(
        p,
        yscale = :log10,
        ylabel = L"\epsilon_u",
        xlabel = "t",
        legend = :bottomleft,
    )

    p
end
function plot_convergence(results, s = :error_u)
    p = plot()

    for result in results
        xs = 5 * result.scales
        plot!(
            p,
            xs,
            getfield.(result.results, s),
            label=string(result.quadrature),
            line = line_style(result.quadrature),
            # marker = marker_style(result.quadrature)
        )
    end

    xs = 5 * results[1].scales
    # plot!(p, 1:results[1].scales[end], x -> 1E-2 * x.^(-1), label=L"\mathcal{O}(x^{-1})", linecolor = :gray)
    # plot!(p, xs, x -> 5E-3 * x.^(-1.5), label=L"\mathcal{O}(x^{-1.5})", linecolor = :gray, linealpha = 0.2, linestyle = :dash)
    # plot!(p, 1:results[1].scales[end], x -> 1E-4 * x.^(-2), label=L"\mathcal{O}(x^{-2})", linecolor = :gray)
    # plot!(p, xs, x -> 5E-4 * x.^(-3), label=L"\mathcal{O}(x^{-3})", linecolor = :gray, linealpha = 0.2, linestyle = :dash)

        plot!(p, xs, x -> 1E-4 * x.^(-2), label=L"\mathcal{O}(x^{-2})", linecolor = :blue, linealpha = 0.2, linestyle = :dash)
        plot!(p, xs, x -> 2E-2 * x.^(-1), label=L"\mathcal{O}(x^{-1})", linecolor = :red, linealpha = 0.2, linestyle = :dash)
        plot!(p, xs, x -> 3E-1 * x.^(-0), label=L"\mathcal{O}(x^{-0})", linecolor = :orange, linealpha = 0.2, linestyle = :dash)

    plot!(
        p,
        yscale = :log10,
        xscale = :log10,
        ylabel = L"\epsilon",
        xlabel = "N",
        legend = :bottomleft,
    )

    p
end
function main2()

    τ = 1.0
    scale = 2
    quadratures = [
        D2Q4(),
        D2Q5(),
        D2Q9(),
        D2Q13(),
        D2Q17(),
        D2Q21(),
        D2Q37(),
    ]

    results = map(quadratures) do q
        couette_velocity_profile(
            q,
            ZeroVelocityInitialCondition(),
            τ,
            scale
        )
    end
    plot_error_locations(results) |> display


    results = map(quadratures) do q
        couette_velocity_profile(
            q,
            ZeroVelocityInitialCondition(),
            τ,
            2scale
        )
    end
    plot_error_locations(results) |> display

    results = map(quadratures) do q
        couette_velocity_profile(
            q,
            ZeroVelocityInitialCondition(),
            τ,
            4scale
        )
    end
    plot_error_locations(results) |> display
end


function main(τ = 1.0, scale = 2)
    quadratures = [
        D2Q4(),
        D2Q5(),
        D2Q9(),
        D2Q13(),
        D2Q17(),
        D2Q21(),
        D2Q37(),
    ]

    results = map(quadratures) do q
        couette_velocity_profile(
            q,
            ZeroVelocityInitialCondition(),
            τ,
            scale
        )
    end

    # scales = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    # scales = [1, 2, 4, 8, 12, 16]
    scales = [1, 2, 4, 8, 16, 32, 64]

    scales = [1, 2, 4, 8, 16]
    scales = [1, 2, 4, 8]

    iteration_strategies = [
        ZeroVelocityInitialCondition(),
        IterativeInitializationMeiEtAl(τ, 1E-7),
        AnalyticalEquilibrium(),
        AnalyticalEquilibriumAndOffEquilibrium(),
    ]

    convergence_results = map(quadratures) do q
        couette_convergence_analysis(
            q,
            iteration_strategies[1],
            τ,
            scales = scales
        )
    end

    convergence_results_iterative = map(quadratures) do q
        couette_convergence_analysis(
            q,
            iteration_strategies[2],
            τ,
            scales = scales
        )
    end
    convergence_results_equilibrium = map(quadratures) do q
        couette_convergence_analysis(
            q,
            iteration_strategies[3],
            τ,
            scales = scales
        )
    end
    convergence_results_offequilibrium = map(quadratures) do q
        couette_convergence_analysis(
            q,
            iteration_strategies[4],
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
    convergence = plot_convergence(main.convergence_results)
    convergence_iterative = plot_convergence(main.convergence_results_iterative)
    convergence_equilibrium = plot_convergence(main.convergence_results_equilibrium)
    convergence_offequilibrium = plot_convergence(main.convergence_results_offequilibrium)

    display(error_location)
    display(error_progression)
    display(convergence)
    display(convergence_iterative)
    display(convergence_equilibrium)
    display(convergence_offequilibrium)


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

end
end

# plot(
#     plot(
#         plot(Examples.Couette.plot_convergence(biggest_result.convergence_results, :error_u), title="Zero"),
#         plot(Examples.Couette.plot_convergence(biggest_result.convergence_results_iterative, :error_u), title = "iterative"),
#         plot(Examples.Couette.plot_convergence(biggest_result.convergence_results_equilibrium, :error_u), title = "equilibrium"),
#         plot(Examples.Couette.plot_convergence(biggest_result.convergence_results_offequilibrium, :error_u), title = "offequilibrium"),
#         legend = nothing,
#     ),

#     plot(
#         # plot(Examples.Couette.plot_convergence(biggest_result.convergence_results, :error_p), title="Zero"),
#         plot(Examples.Couette.plot_convergence(biggest_result.convergence_results_iterative, :error_p), title = "iterative"),
#         plot(Examples.Couette.plot_convergence(biggest_result.convergence_results_equilibrium, :error_p), title = "equilibrium"),
#         plot(Examples.Couette.plot_convergence(biggest_result.convergence_results_offequilibrium, :error_p), title = "offequilibrium"),
#         legend = nothing,
#     ),

#     plot(
#         plot(Examples.Couette.plot_convergence(biggest_result.convergence_results, :error_σ_xy), title="Zero"),
#         plot(Examples.Couette.plot_convergence(biggest_result.convergence_results_iterative, :error_σ_xy), title = "iterative"),
#         plot(Examples.Couette.plot_convergence(biggest_result.convergence_results_equilibrium, :error_σ_xy), title = "equilibrium"),
#         plot(Examples.Couette.plot_convergence(biggest_result.convergence_results_offequilibrium, :error_σ_xy), title = "offequilibrium"),
#         legend = nothing,
#     ),
#     size=(1200, 900)
# )
