module Examples
module TGV_Example

using LatticeBoltzmann, Plots, DataFrames
using LaTeXStrings

using ElectronDisplay
using LaTeXStrings
using JLD2

# using PGFPlots
# using PGFPlotsX
# pgfplots() # x?

import LatticeBoltzmann:
    StopCriteria,
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
    ShowVelocityError,
    D2Q4,
    D2Q5,
    D2Q9,
    D2Q13,
    D2Q17,
    D2Q21,
    D2Q37,
    Quadrature

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
isnothing(x) = false
isnothing(x::Nothing) = true

function tgv_velocity_profile(
    q = D2Q9(),
    initialization_strategy = AnalyticalEquilibrium(),
    τ = 1.0,
    scale = 2,
)
    ν = τ / (2.0 * q.speed_of_sound_squared)
    problem = TaylorGreenVortex(ν, scale, static = true)
    problem = LatticeBoltzmann.TGV(q, τ, scale)

    @show q initialization_strategy

    t_end = 0.20
    plot_every = 100

    if typeof(initialization_strategy) == ZeroVelocityInitialCondition
        @warn "HUH"
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

function tgv_convergence_analysis(
    q = D2Q9(),
    initialization_strategy = AnalyticalEquilibrium(),
    τ = 1.0;
    scales = [1, 2, 4, 8],
)
    ν = τ / (2.0 * q.speed_of_sound_squared)

    t_end = 0.20
    if typeof(initialization_strategy) == ZeroVelocityInitialCondition
        t_end = 3.0
    end
    t_end = 10.0

    results = map(scales) do scale
        # problem = TaylorGreenVortex(ν, scale, static = true)
        problem = LatticeBoltzmann.TGV(q, τ, scale, 16 * scale, 8 * scale)

        Δt = delta_t(problem)
        n_steps = round(Int, t_end / Δt)

        # CompareWithAnalyticalSolution
        process_method =
                TrackHydrodynamicErrors(problem, false, n_steps, StopCriteria(problem))
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

    return (quadrature = q, scales = scales, results = results)

    @show results

    p = (
        plot(
            map(scale -> 8 * scale, scales),
            getfield.(results, :error_u),
            scale = :log10,
            title = "u",
        ),
        plot(
            map(scale -> 8 * scale, scales),
            getfield.(results, :error_σ_xx),
            scale = :log10,
            title = "sigma xx",
        ),
        plot(
            map(scale -> 8 * scale, scales),
            getfield.(results, :error_σ_xy),
            scale = :log10,
            title = "sigma xy",
        ),
        plot(
            map(scale -> 8 * scale, scales),
            getfield.(results, :error_p),
            scale = :log10,
            title = "p",
        ),
    )

    p = plot(
        map(scale -> 8 * scale, scales),
        getfield.(results, :error_u),
        scale = :log10,
        label = "u",
    )
    # plot!(p, map(scale -> 8*scale, scales), getfield.(results, :error_σ_xx), scale=:log10, title="sigma xx")
    plot!(
        p,
        map(scale -> 8 * scale, scales),
        getfield.(results, :error_σ_xy),
        scale = :log10,
        label = "sigma xy",
    )
    plot!(
        p,
        map(scale -> 8 * scale, scales),
        getfield.(results, :error_p),
        scale = :log10,
        label = "p",
    )
    plot!(p, x -> 0.1 * x .^ (-2), scale = :log10, label = "x^-2")
    plot!(p, x -> 0.1 * x .^ (-4), scale = :log10, label = "x^-4")

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
        for y_idx in 1:(problem.NY)
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

    plot!(p, legend = :bottomleft, ylabel = L"\epsilon_u", xlabel = "y", yscale = :log10)

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
            Δt * (1:length(result.processing_method.l2_norms)),
            result.processing_method.l2_norms,
            line = line_style(result.quadrature),
            # marker = marker_style(result.quadrature),
            label = LaTeXString(string(q)),
        )
    end

    plot!(p, yscale = :log10, ylabel = L"\epsilon_u", xlabel = "t", legend = :bottomleft)

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
            label = string(result.quadrature),
            line = line_style(result.quadrature),
            # marker = marker_style(result.quadrature)
        )
    end

    xs = 5 * results[1].scales

    if s == :error_u
        # plot!(p, xs, x -> 1E0 * x.^(-0.5), label=L"\mathcal{O}(x^{-0.5})", linecolor = :gray, linealpha = 0.2, linestyle = :dash)
        # plot!(p, xs, x -> 3E-1 * x.^(-1.5), label=L"\mathcal{O}(x^{-1.5})", linecolor = :gray, linealpha = 0.2, linestyle = :dash)
        # plot!(p, xs, x -> 1E-1 * x.^(-3), label=L"\mathcal{O}(x^{-3})", linecolor = :gray, linealpha = 0.2, linestyle = :dash)

        plot!(
            p,
            xs,
            x -> 1E-1 * x .^ (-2),
            label = L"\mathcal{O}(x^{-2})",
            linecolor = :gray,
            linealpha = 0.2,
            linestyle = :dash,
        )
        plot!(
            p,
            xs,
            x -> 1E-1 * x .^ (-1),
            label = L"\mathcal{O}(x^{-1})",
            linecolor = :gray,
            linealpha = 0.2,
            linestyle = :dash,
        )
        plot!(
            p,
            xs,
            x -> 1E-1 * x .^ (-0),
            label = L"\mathcal{O}(x^{-0})",
            linecolor = :gray,
            linealpha = 0.2,
            linestyle = :dash,
        )
    end

    if s == :error_p
        # plot!(p, xs, x -> 5E-2 * x.^(-2.5), label=L"\mathcal{O}(x^{-2.5})", linecolor = :gray, linealpha = 0.2, linestyle = :dash)
        # plot!(p, xs, x -> 5E-4 * x.^(-5), label=L"\mathcal{O}(x^{-5})", linecolor = :gray, linealpha = 0.2, linestyle = :dash)

        plot!(
            p,
            xs,
            x -> 1E-1 * x .^ (-2),
            label = L"\mathcal{O}(x^{-2})",
            linecolor = :gray,
            linealpha = 0.2,
            linestyle = :dash,
        )
        plot!(
            p,
            xs,
            x -> 1E-3 * x .^ (-4),
            label = L"\mathcal{O}(x^{-4})",
            linecolor = :gray,
            linealpha = 0.2,
            linestyle = :dash,
        )
    end

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
    quadratures = [D2Q4(), D2Q5(), D2Q9(), D2Q13(), D2Q17(), D2Q21(), D2Q37()]

    results = map(quadratures) do q
        tgv_velocity_profile(q, AnalyticalEquilibrium(), τ, scale)
    end
    plot_error_locations(results) |> display

    results = map(quadratures) do q
        tgv_velocity_profile(q, AnalyticalEquilibrium(), τ, 2scale)
    end
    plot_error_locations(results) |> display

    results = map(quadratures) do q
        tgv_velocity_profile(q, AnalyticalEquilibrium(), τ, 4scale)
    end
    plot_error_locations(results) |> display
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

    results = map(quadratures) do q
        tgv_velocity_profile(q, AnalyticalEquilibrium(), τ, scale)
    end

    # scales = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    # scales = [1, 2, 4, 8, 12, 16]
    # scales = [1, 2, 4, 8, 16, 32, 64]

    scales = [1, 2, 4, 8, 16]

    scales = [1 // 2, 1, 2, 4, 8, 16, 32, 64]

    iteration_strategies = [
        IterativeInitializationMeiEtAl(τ, 1E-7),
        AnalyticalEquilibrium(),
        AnalyticalEquilibriumAndOffEquilibrium(),
    ]

    convergence_results = map(quadratures) do q
        tgv_convergence_analysis(q, iteration_strategies[1], τ, scales = scales)
    end

    convergence_results_equilibrium = convergence_results

    convergence_results_iterative = map(quadratures) do q
        tgv_convergence_analysis(q, iteration_strategies[2], τ, scales = scales)
    end
    convergence_results_offequilibrium = map(quadratures) do q
        tgv_convergence_analysis(q, iteration_strategies[3], τ, scales = scales)
    end
    # convergence_results_offequilibrium = map(quadratures) do q
    #     tgv_convergence_analysis(
    #         q,
    #         iteration_strategies[4],
    #         τ,
    #         scales = scales
    #     )
    # end

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
        ),
    )
end

function kruger_analysis()
    q = D2Q9()
    Nx = 31
    Ny = 17

    scales = [1, 2, 4, 8, 16]
    scales = [1, 2, 4]

    Nx = 16
    Ny = 16

    Nx = 31
    Ny = 17
    τs = [0.51, 0.6, 0.8, 1.0]
    τs = [3.0, 2.0, 1.0, 0.8, 0.6, 0.51]
    scales = [1, 2, 4, 8, 16]

    initialization = LatticeBoltzmann.AnalyticalEquilibrium()

    τ_results = map(τs) do τ
        results = map(scales) do scale
            u_0 = sqrt(0.01)
            problem = LatticeBoltzmann.TGV(q, τ, scale, Nx * scale, Ny * scale, u_0 / scale)

            t_end = round(Int, LatticeBoltzmann.decay_time(problem))

            @show t_end

            @time model = LatticeBoltzmann.LatticeBoltzmannModel(
                problem,
                q,
                initialization_strategy = initialization,
                process_method = LatticeBoltzmann.ProcessingMethod(problem, false, t_end),
            )

            @time solution = simulate(model, 1:t_end)
        end

        convergence = map(results) do result
            problem = result.processing_method.problem
            errors = result.processing_method.df[end]

            return (
                τ = problem.τ,
                ν = problem.ν,
                NX = problem.NX,
                NY = problem.NY,
                error_u = errors.error_u,
                error_p = errors.error_p,
                error_σ_xx = errors.error_σ_xx,
                error_σ_xy = errors.error_σ_xy,
            )
        end

        # Show a plot while waiting for the results
        p = plot()
        plot!(
            p,
            getfield.(convergence, :NX),
            getfield.(convergence, :error_u),
            label = L"\epsilon_u",
        )
        plot!(
            p,
            getfield.(convergence, :NX),
            getfield.(convergence, :error_p),
            label = L"\epsilon_p",
        )
        plot!(
            p,
            getfield.(convergence, :NX),
            getfield.(convergence, :error_σ_xx),
            label = L"\epsilon_{\sigma_{xx}}",
        )
        plot!(
            p,
            getfield.(convergence, :NX),
            getfield.(convergence, :error_σ_xy),
            label = L"\epsilon_{\sigma_{xx}}",
        )
        plot!(p, x -> 1E1 * x^(-2), label = L"x^{-2}")
        plot!(p, scale = :log10)
        plot!(p, title = latexstring("\tau = ", τ))
        display(p)

        return results
    end

    convergence = map(τ_results) do τ_result
        map(τ_result) do result
            problem = result.processing_method.problem
            errors = result.processing_method.df[end]

            return (
                τ = problem.τ,
                ν = problem.ν,
                NX = problem.NX,
                NY = problem.NY,
                error_u = errors.error_u,
                error_p = errors.error_p,
                error_σ_xx = errors.error_σ_xx,
                error_σ_xy = errors.error_σ_xy,
            )
        end
    end

    p_u = plot(legend = nothing, scale = :log10, xlabel = L"N", ylabel = L"\epsilon_u")
    p_p = plot(legend = nothing, scale = :log10, xlabel = L"N", ylabel = L"\epsilon_p")
    p_σ_xx = plot(
        legend = nothing,
        scale = :log10,
        xlabel = L"N",
        ylabel = L"\epsilon_{\sigma_{xx}}",
    )
    p_σ_xy = plot(
        legend = :topright,
        scale = :log10,
        xlabel = L"N",
        ylabel = L"\epsilon_{\sigma_{xy}}",
    )
    for τ_convergence in convergence
        label = latexstring("\\tau = ", τ_convergence[1].τ)

        scatter!(
            p_u,
            getfield.(τ_convergence, :NX),
            getfield.(τ_convergence, :error_u),
            label = label,
        )
        scatter!(
            p_p,
            getfield.(τ_convergence, :NX),
            getfield.(τ_convergence, :error_p),
            label = label,
        )
        scatter!(
            p_σ_xx,
            getfield.(τ_convergence, :NX),
            getfield.(τ_convergence, :error_σ_xx),
            label = label,
        )
        scatter!(
            p_σ_xy,
            getfield.(τ_convergence, :NX),
            getfield.(τ_convergence, :error_σ_xy),
            label = label,
        )
    end

    plot!(p_u, x -> 1E0 * x .^ (-2), label = L"\mathcal{O}(x^{-2})", linecolor = :gray)
    plot!(p_p, x -> 1E-1 * x .^ (-2), label = L"\mathcal{O}(x^{-2})", linecolor = :gray)
    plot!(p_σ_xx, x -> 1E0 * x .^ (-2), label = L"\mathcal{O}(x^{-2})", linecolor = :gray)
    plot!(p_σ_xy, x -> 1E0 * x .^ (-2), label = L"\mathcal{O}(x^{-2})", linecolor = :grah)

    plot(p_u, p_p, p_σ_xx, p_σ_xy, size = (900, 600))
    return τ_results

end

function initial_conditions_analysis()
    q = D2Q9()
    scales = [1, 2, 4]

    Nx = 31
    Ny = 17
    τs = [0.51, 0.6, 0.8, 1.0]

    scale = 2

    Nx = 31
    Ny = 17

    τ = 0.8
    t_d = 840

    iteration_strategies = [
        LatticeBoltzmann.AnalyticalVelocityAndStress(),
        LatticeBoltzmann.AnalyticalEquilibriumAndOffEquilibrium(),
    ]
    iteration_strategies = [
        # With this initialization strategy it is assumed that the initial density field
        # (or equivalently the pressure field p) is not available
        LatticeBoltzmann.ConstantDensity(),

        # TODO
        LatticeBoltzmann.AnalyticalVelocityAndStress(),

        # The initial pressure field p is known, which, using the equation fo state,
        # is used to set the initial density
        LatticeBoltzmann.AnalyticalEquilibrium(),

        # Both an initial pressure field and gradient of the velocity is known.
        # The gradient of the velocity is used to initialize the offequilibrium components
        LatticeBoltzmann.AnalyticalEquilibriumAndOffEquilibrium(),

        # We initialize f using an iterative procedure where only the density is conserved
        # it was shown that this procedure gives consistent itnitial conditions for both
        # the equilibrium and off equilibrium components
        LatticeBoltzmann.IterativeInitializationMeiEtAl(τ, 1E-10),
        LatticeBoltzmann.IterativeInitializationMeiEtAl(1.0, 1E-10),
    ]

    init_res = map(iteration_strategies) do initialization
        u_0 = sqrt(0.01)
        u_0 = 0.03
        Nx = 96
        Ny = 72
        τ = 0.8

        # u_0 = 0.06
        # Nx = 48
        # Ny = 36

        # u_0 = 0.01
        # Nx = 32
        # Ny = 32
        # ν = 0.002
        # τ = q.speed_of_sound_squared * ν + 0.5

        problem = LatticeBoltzmann.TGV(q, τ, scale, Nx, Ny, u_0)
        t_end = round(Int, LatticeBoltzmann.decay_time(problem))
        @show t_end

        @time model = LatticeBoltzmann.LatticeBoltzmannModel(
            problem,
            q,
            initialization_strategy = initialization,
            process_method = LatticeBoltzmann.ProcessingMethod(problem, true, t_end),
        )

        @time solution = simulate(model, 1:t_end)
    end

    return init_res
    init_errors = map(init_res) do result
        problem = result.processing_method.problem
        errors = result.processing_method.df

        return (
            error_u = getfield.(errors, :error_u),
            error_p = getfield.(errors, :error_p),
            error_σ_xx = getfield.(errors, :error_σ_xx),
            error_σ_xy = getfield.(errors, :error_σ_xy),
            mass = getfield.(errors, :mass),
            momentum = getfield.(errors, :momentum),
            energy = getfield.(errors, :energy),
        )
    end

    p_u = plot(legend = nothing, scale = :log10, xlabel = L"t", ylabel = L"\epsilon_u")
    p_p = plot(legend = nothing, scale = :log10, xlabel = L"t", ylabel = L"\epsilon_p")
    p_σ_xx = plot(
        legend = nothing,
        scale = :log10,
        xlabel = L"t",
        ylabel = L"\epsilon_{\sigma_{xx}}",
    )
    p_σ_xy = plot(scale = :log10, xlabel = L"t", ylabel = L"\epsilon_{\sigma_{xy}}")

    p_mass = plot(legend = nothing, xlabel = L"t", ylabel = L"Mass")
    p_momentum = plot(legend = nothing, xlabel = L"t", ylabel = L"Momentum")
    p_energy = plot(legend = nothing, xlabel = L"t", ylabel = L"Energy")

    for (errors, init) in zip(
        init_errors,
        [
            "Velocity",
            "Velocity + stress",
            "Velocity + pressure",
            "All",
            "iterative",
            L"iterative, \tau = 1.0",
        ],
    )
        plot!(p_u, errors.error_u, label = init)
        plot!(p_p, errors.error_p, label = init)
        plot!(p_σ_xx, errors.error_σ_xx, label = init)
        plot!(p_σ_xy, errors.error_σ_xy, label = init)

        plot!(p_mass, errors.mass, label = init)
        plot!(p_momentum, errors.momentum, label = init)
        plot!(p_energy, errors.energy, label = init)
    end
    plot(p_u, p_p, p_σ_xx, p_σ_xy) |> display
    plot(p_mass, p_momentum, p_energy) |> display
    return init_res
end

end
end

# plot(
#     plot(
#         plot(Examples.TGV_Example.plot_convergence(result.convergence_results, :error_u), title="Zero"),
#         plot(Examples.TGV_Example.plot_convergence(result.convergence_results_iterative, :error_u), title = "iterative"),
#         plot(Examples.TGV_Example.plot_convergence(result.convergence_results_equilibrium, :error_u), title = "equilibrium"),
#         plot(Examples.TGV_Example.plot_convergence(result.convergence_results_offequilibrium, :error_u), title = "offequilibrium"),
#         legend = nothing,
#     ),

#     plot(
#         plot(Examples.TGV_Example.plot_convergence(result.convergence_results, :error_p), title="Zero"),
#         plot(Examples.TGV_Example.plot_convergence(result.convergence_results_iterative, :error_p), title = "iterative"),
#         plot(Examples.TGV_Example.plot_convergence(result.convergence_results_equilibrium, :error_p), title = "equilibrium"),
#         plot(Examples.TGV_Example.plot_convergence(result.convergence_results_offequilibrium, :error_p), title = "offequilibrium"),
#         legend = nothing,
#     ),

#     plot(
#         plot(Examples.TGV_Example.plot_convergence(result.convergence_results, :error_σ_xx), title="Zero"),
#         plot(Examples.TGV_Example.plot_convergence(result.convergence_results_iterative, :error_σ_xx), title = "iterative"),
#         plot(Examples.TGV_Example.plot_convergence(result.convergence_results_equilibrium, :error_σ_xx), title = "equilibrium"),
#         plot(Examples.TGV_Example.plot_convergence(result.convergence_results_offequilibrium, :error_σ_xx), title = "offequilibrium"),
#         legend = nothing,
#     ),
#     # plot(
#     #     plot(Examples.TGV_Example.plot_convergence(result.convergence_results, :error_σ_xy), title="Zero"),
#     #     plot(Examples.TGV_Example.plot_convergence(result.convergence_results_iterative, :error_σ_xy), title = "iterative"),
#     #     plot(Examples.TGV_Example.plot_convergence(result.convergence_results_equilibrium, :error_σ_xy), title = "equilibrium"),
#     #     plot(Examples.TGV_Example.plot_convergence(result.convergence_results_offequilibrium, :error_σ_xy), title = "offequilibrium"),
#     #     legend = nothing,
#     # ),
#     size=(1200, 900)
# )
