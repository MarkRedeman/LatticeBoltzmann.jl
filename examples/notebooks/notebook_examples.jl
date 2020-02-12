using LatticeBoltzmann, Plots, LaTeXStrings, DataFrames, JLD2
import LatticeBoltzmann: StopCriteria,
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

function shear_wave_convergence_analysis(
    q = D2Q9(),
    initialization_strategy = AnalyticalEquilibrium(),
    τ_lb = 1.0;
    scales = [1, 2, 4, 8],
)
    ν_lb = τ_lb / (2.0 * q.speed_of_sound_squared)
    t_end = 1.0

    results = map(scales) do scale
        problem = DecayingShearFlow(ν_lb, scale, static = true)

        Δt = delta_t(problem)
        n_steps = round(Int, t_end / Δt)

        # CompareWithAnalyticalSolution
        process_method = TrackHydrodynamicErrors(
            problem,
            false,
            n_steps,
            LatticeBoltzmann.NoStoppingCriteria(),
        )

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
end

function plot_snapshots(problem, snapshot_results, snapshots_at, q)
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

    for (time, f_in) in zip(snapshots_at, snapshot_results)
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

            σ_ = LatticeBoltzmann.dimensionless_stress(problem, σ_)
            σ_xx[x_idx, y_idx] = σ_[1, 1]
            σ_xy[x_idx, y_idx] = σ_[1, 2]
        end
        s = (1000, 500)
        y_pos = round(Int, problem.NY / 2)
        domain = x_range[1:Nx]

        x_range, y_range = range(problem)
        velocity = (x, y, t) -> LatticeBoltzmann.velocity(problem, x, y, t)
        σ = (x, y, t) -> LatticeBoltzmann.deviatoric_tensor(q, problem, x, y, t)
        pr = (x, y, t) -> LatticeBoltzmann.pressure(q, problem, x, y, t)

        exact_range = range(0.0, length = 1000, stop = problem.domain_size[1])

        scatter!(
            velocity_profile_x,
            domain,
            u[:, y_pos, 1],
            label = latexstring("t = ", time),
        )
        plot!(
            (x) -> velocity(x, y_range[y_pos], time)[1],
            exact_range,
            label = "",
            linecolor = :gray,
            linealpha = 0.2,
            linestyle = :dash,
        )

        scatter!(
            velocity_profile_y,
            domain,
            u[:, y_pos, 2],
            label = latexstring("t = ", time),
            markershape = :auto,
            markersize = 6,
        )
        plot!(
            (x) -> velocity(x, y_range[y_pos], time)[2],
            exact_range,
            label = "",
            linecolor = :gray,
            linealpha = 0.2,
            linestyle = :dash,
        )

        scatter!(pressure_profile, domain, p[:, y_pos], label = latexstring("t = ", time))
        plot!(
            pressure_profile,
            (x) -> pr(x, y_range[y_pos], time),
            exact_range,
            label = "",
            linecolor = :gray,
            linealpha = 0.2,
            linestyle = :dash,
        )

        plot!(temperature_profile, domain, T[:, y_pos], label = latexstring("t = ", time))

        scatter!(
            sigma_xx_profile,
            domain,
            σ_xx[:, y_pos],
            label = latexstring("t = ", time),
        )
        plot!(
            sigma_xx_profile,
            (x) -> σ(x, y_range[y_pos], time)[1, 1],
            exact_range,
            label = "",
            linecolor = :gray,
            linealpha = 0.2,
            linestyle = :dash,
        )

        scatter!(
            sigma_xy_profile,
            domain,
            σ_xy[:, y_pos],
            label = latexstring("t = ", time),
        )
        plot!(
            sigma_xy_profile,
            (x) -> σ(x, y_range[y_pos], time)[1, 2],
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

function plot_snapshots(problem::CouetteFlow, snapshot_results, snapshots_at, q)
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

    for (time, f_in) in zip(snapshots_at, snapshot_results)
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

        scatter!(
            velocity_profile_x,
            domain,
            u[x_pos, :, 1],
            label = latexstring("t = ", time),
        )
        plot!(
            (y) -> velocity(x_range[x_pos], y, time)[1],
            exact_range,
            label = "",
            linecolor = :gray,
            linealpha = 0.2,
            linestyle = :dash,
        )

        scatter!(
            velocity_profile_y,
            domain,
            u[x_pos, :, 1],
            label = latexstring("t = ", time),
            markershape = :auto,
            markersize = 6,
        )
        plot!(
            velocity_profile_y,
            (y) -> velocity(x_range[x_pos], y, time)[1],
            exact_range,
            label = "",
            linecolor = :gray,
            linealpha = 0.2,
            linestyle = :dash,
        )

        scatter!(pressure_profile, domain, p[x_pos, :], label = latexstring("t = ", time))
        plot!(
            pressure_profile,
            (y) -> pr(x_range[x_pos], y, time),
            exact_range,
            label = "",
            linecolor = :gray,
            linealpha = 0.2,
            linestyle = :dash,
        )

        plot!(temperature_profile, domain, T[x_pos, :], label = latexstring("t = ", time))

        scatter!(
            sigma_xx_profile,
            domain,
            σ_xx[x_pos, :],
            label = latexstring("t = ", time),
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
            label = latexstring("t = ", time),
        )
        # plot!(
        #     sigma_xy_profile,
        #     (y) -> σ(x_range[x_pos], y, time)[1, 2],
        #     exact_range,
        #     label = "",
        #     linecolor = :gray,
        #     linealpha = 0.2,
        #     linestyle = :dash,
        # )
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

function plot_snapshots(problem::PoiseuilleFlow, snapshot_results, snapshots_at, q)
    velocity_profile_x = plot(xlabel = "x", ylabel = L"u_x")
    velocity_profile_y = plot(
        xlabel = L"u_x",
        ylabel = L"y",
        legend = :bottomright,
        title = "Velocity profile at y = pi",
    )
    pressure_profile = plot(xlabel = "y", ylabel = L"p")
    temperature_profile = plot(xlabel = "y", ylabel = L"T")
    sigma_xx_profile = plot(xlabel = "y", ylabel = L"\sigma_{xx}")
    sigma_xy_profile = plot(xlabel = "y", ylabel = L"\sigma_{xy}")

    for (time, f_in) in zip(snapshots_at, snapshot_results)
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

        scatter!(
            velocity_profile_x,
            domain,
            u[x_pos, :, 1],
            label = latexstring("t = ", time),
        )
        plot!(
            velocity_profile_x,
            (y) -> velocity(x_range[x_pos], y, time)[1],
            exact_range,
            label = "",
            linecolor = :gray,
            linealpha = 0.2,
            linestyle = :dash,
        )

        scatter!(
            velocity_profile_y,
            u[x_pos, :, 1],
            domain,
            label = latexstring("t = ", time),
            markershape = :auto,
            markersize = 6,
        )
        plot!(
            velocity_profile_y,
            map((y) -> velocity(x_range[x_pos], y, time)[1], exact_range),
            exact_range,
            label = "",
            linecolor = :gray,
            linealpha = 0.2,
            linestyle = :dash,
        )

        scatter!(pressure_profile, domain, p[x_pos, :], label = latexstring("t = ", time))
        plot!(
            pressure_profile,
            (y) -> pr(x_range[x_pos], y, time),
            exact_range,
            label = "",
            linecolor = :gray,
            linealpha = 0.2,
            linestyle = :dash,
        )

        plot!(temperature_profile, domain, T[x_pos, :], label = latexstring("t = ", time))

        scatter!(
            sigma_xx_profile,
            domain,
            σ_xx[x_pos, :],
            label = latexstring("t = ", time),
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
            label = latexstring("t = ", time),
        )
        # plot!(
        #     sigma_xy_profile,
        #     (y) -> σ(x_range[x_pos], y, time)[1, 2],
        #     exact_range,
        #     label = "",
        #     linecolor = :gray,
        #     linealpha = 0.2,
        #     linestyle = :dash,
        # )
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
