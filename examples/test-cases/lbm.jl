using LatticeBoltzmann, Plots, DataFrames
using LatticeBoltzmann, Plots, DataFrames
using LaTeXStrings
using JLD2
import LatticeBoltzmann:
    StopCriteria, CompareWithAnalyticalSolution, TrackHydrodynamicErrors
import LatticeBoltzmann:
    StopCriteria,
    CompareWithAnalyticalSolution,
    TrackHydrodynamicErrors,
    ZeroVelocityInitialCondition,
    IterativeInitializationMeiEtAl,
    ConstantDensity,
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

function plot_error_locations(results)
    p = plot()

    for result in results
        f = result.f_stream
        problem = result.processing_method.problem
        q = result.quadrature

        analytical_velocity = (x) -> velocity(problem, x, 1.0)[2]
        x_domain = (0.0, problem.domain_size[1])
        y_domain = (0.0, problem.domain_size[2])

        y_idx = 1
        v_y = zeros(problem.NX)
        v_e = zeros(problem.NX)
        v_a = zeros(problem.NX)
        u = zeros(dimension(q))
        x_range, y_range = range(problem)
        for x_idx in 1:(problem.NX)
            ρ = density(q, f[x_idx, y_idx, :])
            velocity!(q, f[x_idx, y_idx, :], ρ, u)
            u = dimensionless_velocity(problem, u)
            v_e[x_idx] = abs(u[2] - analytical_velocity(x_range[x_idx]))
            v_y[x_idx] = u[2]
            v_a[x_idx] = analytical_velocity(x_range[x_idx])
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
            Δt * (1:length(getfield.(result.processing_method.df, :error_u))),
            getfield.(result.processing_method.df, :error_u),
            line = line_style(result.quadrature),
            # marker = marker_style(result.quadrature),
            label = LaTeXString(string(q)),
        )
    end

    plot!(p, yscale = :log10, ylabel = L"\epsilon_u", xlabel = "t", legend = :bottomleft)

    p
end
function plot_convergence(results, s = :error_u, N_0 = 5)
    p = plot()

    for result in results
        xs = N_0 * result.scales
        plot!(
            p,
            xs,
            getfield.(result.results, s),
            label = string(result.quadrature),
            line = line_style(result.quadrature),
            # marker = marker_style(result.quadrature)
        )
    end

    xs = N_0 * results[1].scales
    # plot!(p, 1:results[1].scales[end], x -> 1E-2 * x.^(-1), label=L"\mathcal{O}(x^{-1})", linecolor = :gray)
    # plot!(p, xs, x -> 5E-3 * x.^(-1.5), label=L"\mathcal{O}(x^{-1.5})", linecolor = :gray, linealpha = 0.2, linestyle = :dash)
    # plot!(p, 1:results[1].scales[end], x -> 1E-4 * x.^(-2), label=L"\mathcal{O}(x^{-2})", linecolor = :gray)
    # plot!(p, xs, x -> 5E-4 * x.^(-3), label=L"\mathcal{O}(x^{-3})", linecolor = :gray, linealpha = 0.2, linestyle = :dash)

    if s == :error_u
        plot!(
            p,
            xs,
            x -> 5E-1 * x .^ (-2),
            label = L"\mathcal{O}(x^{-2})",
            linecolor = :blue,
            linealpha = 0.2,
            linestyle = :dash,
        )
        plot!(
            p,
            xs,
            x -> 3E-1 * x .^ (-0),
            label = L"\mathcal{O}(x^{-0})",
            linecolor = :orange,
            linealpha = 0.2,
            linestyle = :dash,
        )
    end

    if s == :error_σ_xx
        plot!(
            p,
            xs,
            x -> 5E-1 * x .^ (-2),
            label = L"\mathcal{O}(x^{-2})",
            linecolor = :blue,
            linealpha = 0.2,
            linestyle = :dash,
        )
        plot!(
            p,
            xs,
            x -> 3E-1 * x .^ (-0),
            label = L"\mathcal{O}(x^{-0})",
            linecolor = :orange,
            linealpha = 0.2,
            linestyle = :dash,
        )
    end

    if s == :error_p
        plot!(
            p,
            xs,
            x -> 1E-2 * x .^ (-2),
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

function plot_convergence_of_results(results, s = :error_u, N_0 = 5)
    p = plot()

    for result in results
        for s in [:error_u, :error_p, :error_σ_xy]
            xs = N_0 * result.scales
            plot!(
                p,
                xs,
                getfield.(result.results, s),
                label = string(result.quadrature),
                line = line_style(result.quadrature),
                # marker = marker_style(result.quadrature)
            )
        end
    end

    xs = N_0 * results[1].scales
    # plot!(p, 1:results[1].scales[end], x -> 1E-2 * x.^(-1), label=L"\mathcal{O}(x^{-1})", linecolor = :gray)
    # plot!(p, xs, x -> 5E-3 * x.^(-1.5), label=L"\mathcal{O}(x^{-1.5})", linecolor = :gray, linealpha = 0.2, linestyle = :dash)
    # plot!(p, 1:results[1].scales[end], x -> 1E-4 * x.^(-2), label=L"\mathcal{O}(x^{-2})", linecolor = :gray)
    # plot!(p, xs, x -> 5E-4 * x.^(-3), label=L"\mathcal{O}(x^{-3})", linecolor = :gray, linealpha = 0.2, linestyle = :dash)

    # if s == :error_u
    plot!(
        p,
        xs,
        x -> 5E-1 * x .^ (-2),
        label = L"\mathcal{O}(x^{-2})",
        linecolor = :blue,
        linealpha = 0.2,
        linestyle = :dash,
    )
    plot!(
        p,
        xs,
        x -> 3E-1 * x .^ (-0),
        label = L"\mathcal{O}(x^{-0})",
        linecolor = :orange,
        linealpha = 0.2,
        linestyle = :dash,
    )
    # end

    # if s == :error_σ_xx
    #     plot!(p, xs, x -> 5E-1 * x.^(-2), label=L"\mathcal{O}(x^{-2})", linecolor = :blue, linealpha = 0.2, linestyle = :dash)
    #     plot!(p, xs, x -> 3E-1 * x.^(-0), label=L"\mathcal{O}(x^{-0})", linecolor = :orange, linealpha = 0.2, linestyle = :dash)
    # end

    # if s == :error_p
    #     plot!(p, xs, x -> 1E-2 * x.^(-2), label=L"\mathcal{O}(x^{-2})", linecolor = :gray, linealpha = 0.2, linestyle = :dash)
    #     plot!(p, xs, x -> 1E-3 * x.^(-4), label=L"\mathcal{O}(x^{-4})", linecolor = :gray, linealpha = 0.2, linestyle = :dash)
    # end

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

function plot_convergence_of_result(result, s = :error_u, N_0 = 5, q = D2Q9())
    p = plot()

    xs = N_0 * result.scales
    plot!(
        p,
        xs,
        getfield.(result.results, :error_u),
        line = line_style(q),
        label = L"\epsilon_u",
        # marker = marker_style(result.quadrature)
    )
    plot!(
        p,
        xs,
        getfield.(result.results, :error_p),
        line = line_style(q),
        label = L"\epsilon_p",
        # marker = marker_style(result.quadrature)
    )
    plot!(
        p,
        xs,
        getfield.(result.results, :error_σ_xy),
        line = line_style(q),
        label = L"\epsilon_{\sigma_{xy}}",
        # marker = marker_style(result.quadrature)
    )

    xs = N_0 * result.scales
    # plot!(p, 1:results[1].scales[end], x -> 1E-2 * x.^(-1), label=L"\mathcal{O}(x^{-1})", linecolor = :gray)
    # plot!(p, xs, x -> 5E-3 * x.^(-1.5), label=L"\mathcal{O}(x^{-1.5})", linecolor = :gray, linealpha = 0.2, linestyle = :dash)
    # plot!(p, 1:results[1].scales[end], x -> 1E-4 * x.^(-2), label=L"\mathcal{O}(x^{-2})", linecolor = :gray)
    # plot!(p, xs, x -> 5E-4 * x.^(-3), label=L"\mathcal{O}(x^{-3})", linecolor = :gray, linealpha = 0.2, linestyle = :dash)

    # if s == :error_u
    plot!(
        p,
        xs,
        x -> 5E-1 * x .^ (-2),
        label = L"\mathcal{O}(x^{-2})",
        linecolor = :blue,
        linealpha = 0.2,
        linestyle = :dash,
    )
    # end

    # if s == :error_σ_xx
    #     plot!(p, xs, x -> 5E-1 * x.^(-2), label=L"\mathcal{O}(x^{-2})", linecolor = :blue, linealpha = 0.2, linestyle = :dash)
    #     plot!(p, xs, x -> 3E-1 * x.^(-0), label=L"\mathcal{O}(x^{-0})", linecolor = :orange, linealpha = 0.2, linestyle = :dash)
    # end

    # if s == :error_p
    #     plot!(p, xs, x -> 1E-2 * x.^(-2), label=L"\mathcal{O}(x^{-2})", linecolor = :gray, linealpha = 0.2, linestyle = :dash)
    #     plot!(p, xs, x -> 1E-3 * x.^(-4), label=L"\mathcal{O}(x^{-4})", linecolor = :gray, linealpha = 0.2, linestyle = :dash)
    # end

    plot!(
        p,
        yscale = :log10,
        xscale = :log10,
        ylabel = L"\epsilon",
        xlabel = "N",
        title = "Convergence analysis of a steady-state shear wave",
        legend = :bottomleft,
    )

    p
end
