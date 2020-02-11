using LatticeBoltzmann, Plots, DataFrames
using LaTeXStrings

module Examples
module Couette

using LatticeBoltzmann, Plots, DataFrames
using ElectronDisplay
using LaTeXStrings
using JLD2

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
    ShowVelocityError

# line_style(q::Quadrature) = (:line, :dot)

line_style(q::D2Q4) = (:line)
line_style(q::D2Q5) = (:line, :dash, 1.0, 1.0)
line_style(q::D2Q9) = (:line)
line_style(q::D2Q13) = (:line)
line_style(q::D2Q17) = (:line)
line_style(q::D2Q21) = (:line)
line_style(q::D2Q37) = (:line, :dash, 1.0, 1.0)
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
    scale = 2,
)
    ν = (τ - 0.5) / q.speed_of_sound_squared
    ν = 1.0
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
    t_end = 1.0

    simulate(
        problem,
        q,
        should_process = false,
        initialization_strategy = initialization_strategy,
        t_end = t_end,
        collision_model = SRT,
        # collision_model = MRT,
        process_method = ShowVelocityError(
            problem,
            plot_every,
            Float64[],
            LatticeBoltzmann.VelocityConvergenceStoppingCriteria(1E-7, problem),
        ),
    )
end

function couette_convergence_analysis(
    q = D2Q9(),
    initialization_strategy = ZeroVelocityInitialCondition(),
    τ = 1.0;
    scales = [1, 2, 4, 8],
)
    ν = τ / (2.0 * q.speed_of_sound_squared)

    t_end = 0.20
    if typeof(initialization_strategy) == ZeroVelocityInitialCondition
        t_end = 3.0
    end
    t_end = 1.0

    results = map(scales) do scale
        problem = CouetteFlow(ν, scale)

        Δt = delta_t(problem)
        n_steps = round(Int, t_end / Δt)

        # CompareWithAnalyticalSolution
        process_method = TrackHydrodynamicErrors(
            problem,
            false,
            n_steps,
            # LatticeBoltzmann.NoStoppingCriteria()
            LatticeBoltzmann.VelocityConvergenceStoppingCriteria(1E-7, problem),
            # LatticeBoltzmann.MeanVelocityStoppingCriteria(0.0, 1e-7, problem)
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

    plot!(p, yscale = :log10, ylabel = L"\epsilon_u", xlabel = "t", legend = :bottomright)

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
    # plot!(p, 1:results[1].scales[end], x -> 1E-2 * x.^(-1), label=L"\mathcal{O}(x^{-1})", linecolor = :gray)
    # plot!(p, xs, x -> 5E-3 * x.^(-1.5), label=L"\mathcal{O}(x^{-1.5})", linecolor = :gray, linealpha = 0.2, linestyle = :dash)
    # plot!(p, 1:results[1].scales[end], x -> 1E-4 * x.^(-2), label=L"\mathcal{O}(x^{-2})", linecolor = :gray)
    # plot!(p, xs, x -> 5E-4 * x.^(-3), label=L"\mathcal{O}(x^{-3})", linecolor = :gray, linealpha = 0.2, linestyle = :dash)

    plot!(
        p,
        xs,
        x -> 1E-4 * x .^ (-2),
        label = L"\mathcal{O}(x^{-2})",
        linestyle = :dash,
        linecolor = :gray,
    )
    plot!(
        p,
        xs,
        x -> 2E-2 * x .^ (-1),
        label = L"\mathcal{O}(x^{-1})",
        linestyle = :dash,
        linecolor = :gray,
    )
    # plot!(p, xs, x -> 3E-1 * x.^(-0), label=L"\mathcal{O}(x^{-0})", linecolor = :orange, linealpha = 0.2, linestyle = :dash)

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
        couette_velocity_profile(q, ZeroVelocityInitialCondition(), τ, scale)
    end
    plot_error_locations(results) |> display

    results = map(quadratures) do q
        couette_velocity_profile(q, ZeroVelocityInitialCondition(), τ, 2scale)
    end
    plot_error_locations(results) |> display

    results = map(quadratures) do q
        couette_velocity_profile(q, ZeroVelocityInitialCondition(), τ, 4scale)
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
        couette_velocity_profile(
            q,
            # ZeroVelocityInitialCondition(),
            AnalyticalEquilibrium(),
            τ,
            scale,
        )
    end
    plot_error_progresion(results) |> display

    scales = [1, 2, 4]

    iteration_strategies = [
        # AnalyticalEquilibrium(),
        ZeroVelocityInitialCondition(),
        IterativeInitializationMeiEtAl(τ, 1E-7),
        AnalyticalEquilibrium(),
        AnalyticalEquilibriumAndOffEquilibrium(),
    ]

    τ = 0.8
    convergence_results = map(quadratures) do q
        couette_convergence_analysis(q, iteration_strategies[1], τ, scales = scales)
    end

    convergence_results_iterative = map(quadratures) do q
        couette_convergence_analysis(q, iteration_strategies[2], τ, scales = scales)
    end
    convergence_results_equilibrium = map(quadratures) do q
        couette_convergence_analysis(q, iteration_strategies[3], τ, scales = scales)
    end
    convergence_results_offequilibrium = map(quadratures) do q
        couette_convergence_analysis(q, iteration_strategies[4], τ, scales = scales)
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
    # error_location = plot_error_locations(main.results)
    # error_progression = plot_error_progresion(main.results)
    convergence = plot_convergence(main.convergence_results)
    convergence_iterative = plot_convergence(main.convergence_results_iterative)
    convergence_equilibrium = plot_convergence(main.convergence_results_equilibrium)
    convergence_offequilibrium = plot_convergence(main.convergence_results_offequilibrium)

    # display(error_location)
    # display(error_progression)
    display(convergence)
    # display(convergence_iterative)
    # display(convergence_equilibrium)
    # display(convergence_offequilibrium)

    return (
        results = main.results,
        convergence_results = main.convergence_results,
        convergence_results_iterative = main.convergence_results_iterative,
        convergence_results_equilibrium = main.convergence_results_equilibrium,
        convergence_results_offequilibrium = main.convergence_results_offequilibrium,
        plots = (
            # error_location = error_location,
            # error_progression = error_progression,
            convergence = convergence,
        ),
    )
end

function show_velocity_profile_x(scale = 2, snapshots_at = [0, 1, 2, 100, 200, 1000])
    velocity_profiles = map([D2Q9(), D2Q13(), D2Q17(), D2Q21(), D2Q37()]) do q
        result = compute_snapshots(
            q,
            scale = scale,
            snapshots_at = snapshots_at,
            initialization_strategy = AnalyticalEquilibrium(),
        )

        p = result.velocity_profile_x
        plot!(p, title = string(q))

        return result.velocity_profile_x
    end

    return velocity_profiles

    plot(velocity_profiles..., legend = nothing, size = (900, 600))
end

function compute_snapshots(
    q = D2Q9();
    scale = 1,
    snapshots_at = [0, 1, 2, 100, 200, 1000],
    initialization_strategy = AnalyticalEquilibrium(),
)
    ν = 1.0 / (2 * q.speed_of_sound_squared)

    problem = CouetteFlow(ν, scale)
    Δt = delta_t(problem)

    snapshot_results = map(snapshots_at) do t_end
        simulation = LatticeBoltzmann.LatticeBoltzmannModel(
            problem,
            q,
            initialization_strategy = initialization_strategy,
            process_method = LatticeBoltzmann.ProcessingMethod(problem, true, 1),
        )

        if t_end == 0
            return simulation
        end

        simulate(simulation, 1:t_end)
    end

    plot_snapshots(problem, snapshot_results, Δt .* snapshots_at, q)
end

function plot_snapshots(problem, snapshot_results, snapshots, q)
    velocity_profile_x = plot(xlabel = "x", ylabel = L"u_x")
    velocity_profile_y = plot(
        xlabel = "x",
        ylabel = L"u_y",
        legend = :bottomright,
        title = string("Velocity profile at ", latexstring("y = \\pi")),
    )
    pressure_profile = plot(xlabel = "x", ylabel = L"p")
    temperature_profile = plot(xlabel = "x", ylabel = L"T")
    sigma_xx_profile = plot(xlabel = "x", ylabel = L"\sigma_{xx}")
    sigma_xy_profile = plot(xlabel = "x", ylabel = L"\sigma_{xy}")

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
        @inbounds for x_idx in 1:Nx, y_idx in 1:Ny
            @inbounds for f_idx in 1:size(f_in, 3)
                f[f_idx] = f_in[x_idx, y_idx, f_idx]
            end
            x = x_range[x_idx]
            y = y_range[y_idx]

            ρ_ = density(q, f)
            velocity!(q, f, ρ_, u_)
            T_ = temperature(q, f, ρ_, u_)
            p_ = pressure(q, f, ρ_, u_)

            τ = q.speed_of_sound_squared * LatticeBoltzmann.lattice_viscosity(problem)
            σ_ = LatticeBoltzmann.deviatoric_tensor(q, τ, f, ρ_, u_)

            ρ_ = LatticeBoltzmann.dimensionless_density(problem, ρ_)
            u_ = LatticeBoltzmann.dimensionless_velocity(problem, u_)
            T_ = LatticeBoltzmann.dimensionless_temperature(q, problem, T_)
            p_ = LatticeBoltzmann.dimensionless_pressure(q, problem, p_)
            ρ[x_idx, y_idx] = ρ_
            u[x_idx, y_idx, :] = u_
            p[x_idx, y_idx] = p_
            T[x_idx, y_idx] = T_

            u[x_idx, y_idx, :] .-= LatticeBoltzmann.velocity(problem, x, y, time)
            σ_ = LatticeBoltzmann.dimensionless_stress(problem, σ_)
            σ_xx[x_idx, y_idx] = σ_[1, 1]
            σ_xy[x_idx, y_idx] = σ_[1, 2]
        end
        s = (1000, 500)
        x_pos = max(round(Int, problem.NX / 2), 1)
        domain = y_range[1:Ny]

        x_range, y_range = range(problem)
        velocity = (x, y, t) -> LatticeBoltzmann.velocity(problem, x, y, t)
        σ = (x, y, t) -> LatticeBoltzmann.deviatoric_tensor(q, problem, x, y, t)
        pr = (x, y, t) -> LatticeBoltzmann.pressure(q, problem, x, y, t)

        exact_range = range(0.0, length = 1000, stop = problem.domain_size[1])

        @show x_pos problem.NX, problem.NY

        timestep = round(Int, time / delta_t(problem))
        # scatter!(velocity_profile_x, domain, u[x_pos, :, 1], label = latexstring("t = ", time))
        # plot!((y) -> velocity(x_range[x_pos], y, time)[1], exact_range, label = "", linecolor = :gray, linealpha = 0.2, linestyle = :dash)

        marker = nothing

        if timestep == 0
            marker = :circle
        elseif timestep == 1
            marker = :cross
        elseif timestep == 2
            marker = :xcross
        end
        scatter!(
            velocity_profile_x,
            domain,
            abs.(u[x_pos, :, 1]),
            label = latexstring("t = ", timestep),
            markershape = marker,
            markercolor = :gray,
            markersize = 8,
        )

        scatter!(
            velocity_profile_y,
            domain,
            u[x_pos, :, 2],
            label = latexstring("t = ", timestep),
            markershape = :auto,
            markersize = 6,
        )
        plot!(
            (y) -> velocity(x_range[x_pos], y, time)[2],
            exact_range,
            label = "",
            linecolor = :gray,
            linealpha = 0.2,
            linestyle = :dash,
        )

        scatter!(
            pressure_profile,
            domain,
            p[x_pos, :],
            label = latexstring("t = ", timestep),
        )
        plot!(
            pressure_profile,
            (y) -> pr(x_range[x_pos], y, time),
            exact_range,
            label = "",
            linecolor = :gray,
            linealpha = 0.2,
            linestyle = :dash,
        )

        plot!(
            temperature_profile,
            domain,
            T[x_pos, :],
            label = latexstring("t = ", timestep),
        )

        scatter!(
            sigma_xx_profile,
            domain,
            σ_xx[x_pos, :],
            label = latexstring("t = ", timestep),
        )
        plot!(
            sigma_xx_profile,
            (y) -> σ(x_range[x_pos], y, time)[1, 1],
            exact_range,
            label = "",
            linecolor = :gray,
            linealpha = 0.2,
            linestyle = :dash,
        )

        scatter!(
            sigma_xy_profile,
            domain,
            σ_xy[x_pos, :],
            label = latexstring("t = ", timestep),
        )
        plot!(
            sigma_xy_profile,
            (y) -> σ(x_range[x_pos], y, time)[1, 2],
            exact_range,
            label = "",
            linecolor = :gray,
            linealpha = 0.2,
            linestyle = :dash,
        )
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

function hoi(scale, ν = 1.0)
    quadratures = [D2Q4(), D2Q5(), D2Q9(), D2Q13(), D2Q17(), D2Q21(), D2Q37()]

    results = map(quadratures) do q
        τ = 0.5 + q.speed_of_sound_squared * ν
        @show τ
        couette_velocity_profile(
            q,
            ZeroVelocityInitialCondition(),
            # AnalyticalEquilibrium(),
            τ,
            scale,
        )
    end
    plot_error_progresion(results)
end

function hoi_snapshot(q = D2Q9(), initialization_strategy = nothing, τ = 1.0, scale = 2)
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
        collision_model = SRT,
        # collision_model = MRT,
        process_method = LatticeBoltzmann.TakeSnapshots(problem, 100),
    )

    problem = CouetteFlow(1.0 / 6.0, 8)
    ν = LatticeBoltzmann.viscosity(problem)
    Δt = delta_t(problem)

    # ν t / h^2 = 1/2
    # Δt * t_i = t = h^2 / (2ν)
    # t_i = t = h^2 / (2ν * Δt)
    # t_i = t = 0.5 * 1.0 / (2ν * Δt)

    problem = CouetteFlow(1.0 / 6.0, 8)

    problem = CouetteFlow(1.0 / 6.0, 32)
    problem = PoiseuilleFlow(1.0 / 6.0, 8)
    ν = LatticeBoltzmann.viscosity(problem)
    Δt = delta_t(problem)
    snapshot_at =
        round.(Int, [
            # 0.0005 * 1.0 / (ν * Δt),
            0.005 * 1.0 / (ν * Δt),
            0.05 * 1.0 / (ν * Δt),
            # 0.5 * 1.0 / (ν * Δt),
            5.0 * 1.0 / (ν * Δt),
        ])

    simulate(
        problem,
        D2Q9(),
        t_end = 5.0 / LatticeBoltzmann.viscosity(problem),
        process_method = LatticeBoltzmann.TakeSnapshots(problem, snapshot_at),
        initialization_strategy = LatticeBoltzmann.ZeroVelocityInitialCondition(),
    )

    simulate(
        problem,
        D2Q9(),
        t_end = 1.0,
        process_method = LatticeBoltzmann.TakeSnapshots(
            CouetteFlow(1.0 / 6.0, 8),
            [0, 1, 5, 10, 50, 100, 200, 500, 1000, 5000, 10000],
        ),
        initialization_strategy = LatticeBoltzmann.ZeroVelocityInitialCondition(),
    )

end

function truncation_versus_numerical(xs)
    q = D2Q13()
    τ = 0.8
    ν = τ / (2.0 * q.speed_of_sound_squared)
    scale = 2

    NX = 1
    NY = 5 * scale
    domain_size = (1.0, 1.0)
    map(xs) do u_max
        problem = CouetteFlow(1.0, u_max / scale, ν, NX, NY, domain_size)

        t_end = 1.0
        Δt = delta_t(problem)
        @show Δt
        n_steps = round(Int, t_end / Δt)

        # CompareWithAnalyticalSolution
        process_method = TrackHydrodynamicErrors(
            problem,
            false,
            n_steps,
            LatticeBoltzmann.VelocityConvergenceStoppingCriteria(1E-7, problem),
        )
        @show problem

        @time res = simulate(
            problem,
            q,
            process_method = process_method,
            initialization_strategy = ZeroVelocityInitialCondition(),
            #initialization_strategy = AnalyticalEquilibrium(),
            t_end = t_end,
        )

        return res.processing_method.df[end].error_u
    end
end

end
end
