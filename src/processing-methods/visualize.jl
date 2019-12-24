function visualize(problem::FluidFlowProblem, quadrature::Quadrature, f_in, time, stats)
    q = quadrature
    # Density
    ρ = lbm.density(quadrature, f_in)
    ρ = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))
    p = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))
    T = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))
    σ_xx = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))
    σ_xy = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))

    Nx = size(f_in, 1)
    Ny = size(f_in, 2)

    # Momentum
    j = lbm.momentum(quadrature, f_in)
    E = lbm.total_energy(quadrature, f_in)
    E_k = lbm.kinetic_energy(quadrature, f_in, ρ, j ./ ρ)
    ϵ = 1.0 #lbm.internal_energy(quadrature, f_in, ρ, j ./ ρ)

    density_field = fill(0.0, Nx, Ny)
    pressure_field = fill(0.0, Nx, Ny)
    velocity_field = fill(0.0, Nx, Ny, lbm.dimension(quadrature))
    σ_xx_field = fill(0.0, Nx, Ny)
    σ_xy_field = fill(0.0, Nx, Ny)

    x_range, y_range = range(problem)

    f = Array{Float64}(undef, size(f_in, 3))
    u = zeros(dimension(q))
    @inbounds for x_idx = 1:Nx, y_idx = 1:Ny
        @inbounds for f_idx = 1:size(f_in, 3)
            f[f_idx] = f_in[x_idx, y_idx, f_idx]
        end
        x = x_range[x_idx]
        y = y_range[y_idx]

        density_field[x_idx, y_idx] = lbm.density(quadrature, problem, x, y, time)
        pressure_field[x_idx, y_idx] = lbm.pressure(quadrature, problem, x, y, time)
        velocity_field[x_idx, y_idx, :] = lbm.velocity(problem, x, y, time)
        σ_exact = deviatoric_tensor(q, problem, x, y, time)
        σ_xx_field[x_idx, y_idx] = σ_exact[1, 1]
        σ_xy_field[x_idx, y_idx] = σ_exact[1, 2]

        ρ[x_idx, y_idx] = density(q, f)
        p[x_idx, y_idx] =
            pressure(q, f, ρ[x_idx, y_idx], j[x_idx, y_idx, :] ./ ρ[x_idx, y_idx])
        T[x_idx, y_idx] =
            temperature(q, f, ρ[x_idx, y_idx], j[x_idx, y_idx, :] ./ ρ[x_idx, y_idx])


        τ = q.speed_of_sound_squared * lbm.lattice_viscosity(problem)
        σ = deviatoric_tensor(q, τ, f, ρ[x_idx, y_idx], j[x_idx, y_idx, :] ./ ρ[x_idx, y_idx])
        σ = dimensionless_stress(problem, σ)
        σ_xx[x_idx, y_idx] = σ[1, 1]
        σ_xy[x_idx, y_idx] = σ[1, 2]

        ρ[x_idx, y_idx] = dimensionless_density(problem, ρ[x_idx, y_idx])
        p[x_idx, y_idx] = dimensionless_pressure(q, problem, p[x_idx, y_idx])
        T[x_idx, y_idx] = dimensionless_temperature(q, problem, T[x_idx, y_idx])
        j[x_idx, y_idx, :] =
            dimensionless_velocity(problem, j[x_idx, y_idx, :] ./ ρ[x_idx, y_idx])
    end

    s = (1000, 500)

    domain = (2:(problem.NY-1)) ./ (problem.NY - 2)

    if (typeof(problem) != DecayingShearFlow)
        x_pos = round(Int, problem.NX / 2)
        domain = y_range[1:Ny]

        velocity_profile_x = plot(
            domain,
            j[x_pos, 1:(problem.NY), 1],
            label = "solution",
            title = "u_x",
            legend = nothing,
        )
        plot!(
            velocity_profile_x,
            domain,
            velocity_field[x_pos, 1:(problem.NY), 1],
            label = "exact",
        )

        # @show j[x_pos, 1:(problem.NY), 1] velocity_field[x_pos, 1:(problem.NY), 1]
        # velocity_profile_x = plot(domain, j[x_pos, 1:(problem.NY), 1] - velocity_field[x_pos, 1:(problem.NY), 1], label="solution", title="u_x", legend=nothing)
        # scatter!(velocity_profile_x, domain, j[x_pos, 1:(problem.NY), 1] - velocity_field[x_pos, 1:(problem.NY), 1], label="solution", title="u_x", legend=nothing)

        velocity_profile_y = plot(
            j[x_pos, 1:(problem.NY), 2],
            domain,
            label = "solution",
            title = "u_y",
            legend = nothing,
        )
        plot!(
            velocity_profile_y,
            velocity_field[x_pos, 1:(problem.NY), 2],
            domain,
            label = "exact",
        )

        pressure_profile = plot(
            domain,
            p[x_pos, 1:(problem.NY)],
            label = "solution",
            title = "p",
            legend = nothing,
        )
        plot!(
            pressure_profile,
            domain,
            pressure_field[x_pos, 1:(problem.NY)],
            label = "exact",
        )

        temperature_profile = plot(
            domain,
            T[x_pos, 1:(problem.NY), 1],
            label = "solution",
            title = "T",
            legend = nothing,
        )
        # plot!(temperature_profile, domain, temperature_field[x_pos, 1:(problem.NY), 1], label="exact")

        sigma_xx_profile = plot(
            domain,
            σ_xx[x_pos, 1:(problem.NY)],
            label = "solution",
            title = "sigma_xx",
            legend = nothing,
        )
        plot!(
            sigma_xx_profile,
            domain,
            σ_xx_field[x_pos, 1:(problem.NY)],
            label = "exact",
        )

        sigma_xy_profile = plot(
            domain,
            σ_xy[x_pos, 1:(problem.NY)],
            label = "solution",
            title = "sigma_xy",
            legend = nothing,
        )
        plot!(
            sigma_xy_profile,
            domain,
            σ_xy_field[x_pos, 1:(problem.NY)],
            label = "exact",
        )
    else
        y_pos = round(Int, problem.NY / 2)
        domain = x_range[1:Nx]

        velocity_profile_x =
            plot(domain, j[:, y_pos, 1] ./ ρ[:, y_pos], label = "solution", title = "u_x")
        plot!(velocity_profile_x, domain, velocity_field[:, y_pos, 1], label = "exact")

        velocity_profile_y =
            plot(j[:, y_pos, 2] ./ ρ[:, y_pos], domain, label = "solution", title = "u_y")
        plot!(velocity_profile_y, velocity_field[:, y_pos, 2], domain, label = "exact")
        # velocity_profile_y = plot(j[:, y_pos, 2] ./ ρ[:, y_pos] - velocity_field[:, y_pos, 2], domain, label="solution", title="u_y")

        pressure_profile =
            plot(domain, p[:, y_pos], label = "solution", title = "p", legend = nothing)
        plot!(pressure_profile, domain, pressure_field[:, y_pos], label = "exact")

        temperature_profile =
            plot(domain, T[:, y_pos], label = "solution", title = "T", legend = nothing)

        sigma_xx_profile = plot(
            domain,
            σ_xx[:, y_pos],
            label = "solution",
            title = "sigma_xx",
            legend = nothing,
        )
        plot!(
            sigma_xx_profile,
            domain,
            σ_xx_field[:, y_pos],
            label = "exact",
        )

        sigma_xy_profile = plot(
            domain,
            σ_xy[:, y_pos],
            label = "solution",
            title = "sigma_xy",
            legend = nothing,
        )
        plot!(
            sigma_xy_profile,
            domain,
            σ_xy_field[:, y_pos],
            label = "exact",
        )
    end


    # kinetic_energy_profile = plot(getfield.(stats, :kinetic_energy), legend=false, title="Kinetic energy")

    plot(
        contour(ρ, title = "Density", fill = true, cbar = true),
        contour(p, title = "Pressure", fill = true),
        contour(T, title = "Temperature", fill = true),
        # contour(pressure_field', title="pressure analytical", fill=true),
        plot!(streamline(j), title = "Computed"),
        plot!(streamline(velocity_field), title = "exact"),
        # heatmap(j[:, :, 1]', fill=true),
        # heatmap(j[:, 2:(problem.NY-1), 2]', fill=true),
        # heatmap(velocity_field[:, :, 1]', fill=true),
        # heatmap(velocity_field[:, :, 2]', fill=true),
        # streamline(velocity_field .- ρ .* j),
        plot(getfield.(stats, :error_u), legend = false, title = "U_e"),
        plot(getfield.(stats, :error_p), legend = false, title = "P_e"),
        plot(getfield.(stats, :error_σ_xx), legend = false, title = "sigma_xx_e"),
        plot(getfield.(stats, :error_σ_yx), legend = false, title = "sigma_yx_e"),
        plot(getfield.(stats, :error_σ_xy), legend = false, title = "sigma_xy_e"),
        plot(getfield.(stats, :error_σ_yy), legend = false, title = "sigma_yy_e"),
        # plot(getfield.(stats, :error_u), legend=false, title="U_e"),
        # plot(getfield.(stats, :error_p), legend=false, title="P_e"),
        # kinetic_energy_profile,
        pressure_profile,
        temperature_profile,
        velocity_profile_x,
        velocity_profile_y,
        sigma_xx_profile,
        sigma_xy_profile,
        size = (1000, 600),
    )
    gui()
end

function streamline(
    j;
    amount_of_arrows = 5,
    step = round(Int, size(j, 1) / amount_of_arrows),
)
    s = (1000, 500)
    velocity_field = contour(
        (j[:, :, 1] .^ 2 .+ j[:, :, 2] .^ 2)',
        # cbar = true,
        fill = true,
        title = "Momentum",
    )
    N = size(j, 1)
    X = [
        i
        for i in range(2, size(j, 1), step = step), j in range(1, size(j, 2), step = step)
    ]
    Y = [
        j
        for i in range(2, size(j, 1), step = step), j in range(1, size(j, 2), step = step)
    ]
    # @show "process: ", u, v,

    quiver!(
        velocity_field,
        X,
        Y,
        quiver = (x, y) -> (
            0.1 * amount_of_arrows * j[Int(x), Int(y), 1] /
            sqrt(j[Int(x), Int(y), 1]^2 + j[Int(x), Int(y), 2]^2),
            0.1 * amount_of_arrows * j[Int(x), Int(y), 2] /
            sqrt(j[Int(x), Int(y), 1]^2 + j[Int(x), Int(y), 2]^2),
        ),
        color = "white",
    )

    return velocity_field
end
