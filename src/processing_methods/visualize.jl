function visualize(problem::FluidFlowProblem, quadrature::Quadrature, f_in, time, stats)
    return
    q = quadrature

    # Pre-allocate Macroscopic Variables
    ρ = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))
    u = Array{Float64}(undef, size(f_in, 1), size(f_in, 2), dimension(q))
    p = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))
    T = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))
    σ_xx = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))
    σ_xy = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))

    Nx = size(f_in, 1)
    Ny = size(f_in, 2)

    density_field = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))
    velocity_field = Array{Float64}(undef, size(f_in, 1), size(f_in, 2), dimension(q))
    pressure_field = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))
    σ_xx_field = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))
    σ_xy_field = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))

    x_range, y_range = range(problem)

    f = Array{Float64}(undef, size(f_in, 3))
    u_ = zeros(dimension(q))
    @inbounds for x_idx in 1:Nx, y_idx in 1:Ny
        @inbounds for f_idx in 1:size(f_in, 3)
            f[f_idx] = f_in[x_idx, y_idx, f_idx]
        end

        x = x_range[x_idx]
        y = y_range[y_idx]

        density_field[x_idx, y_idx] =
            LatticeBoltzmann.density(quadrature, problem, x, y, time)
        pressure_field[x_idx, y_idx] =
            LatticeBoltzmann.pressure(quadrature, problem, x, y, time)
        velocity_field[x_idx, y_idx, :] = LatticeBoltzmann.velocity(problem, x, y, time)
        σ_exact = deviatoric_tensor(q, problem, x, y, time)
        σ_xx_field[x_idx, y_idx] = σ_exact[1, 1]
        σ_xy_field[x_idx, y_idx] = σ_exact[1, 2]

        ρ_ = density(q, f)
        velocity!(q, f, ρ_, u_)
        T_ = temperature(q, f, ρ_, u_)
        p_ = pressure(q, f, ρ_, u_)

        τ = q.speed_of_sound_squared * LatticeBoltzmann.lattice_viscosity(problem)
        σ_ = deviatoric_tensor(q, τ, f, ρ_, u_)

        ρ_ = dimensionless_density(problem, ρ_)
        u_ = dimensionless_velocity(problem, u_)
        T_ = dimensionless_temperature(q, problem, T_)
        p_ = dimensionless_pressure(q, problem, p_)
        ρ[x_idx, y_idx] = ρ_
        u[x_idx, y_idx, :] = u_
        p[x_idx, y_idx] = p_
        T[x_idx, y_idx] = T_

        σ_ = dimensionless_stress(problem, σ_)
        σ_xx[x_idx, y_idx] = σ_[1, 1]
        σ_xy[x_idx, y_idx] = σ_[1, 2]
    end

    s = (1000, 500)

    domain = (2:(problem.NY - 1)) ./ (problem.NY - 2)

    if (typeof(problem) != DecayingShearFlow)
        x_pos = round(Int, problem.NX / 2)
        domain = y_range[1:Ny]

        velocity_profile_x = plot(domain, u[x_pos, :, 1], label = "solution", title = "u_x")
        plot!(velocity_profile_x, domain, velocity_field[x_pos, :, 1], label = "exact")

        velocity_profile_y =
            plot(u[x_pos, 1:(problem.NY), 2], domain, label = "solution", title = "u_y")
        plot!(
            velocity_profile_y,
            velocity_field[x_pos, 1:(problem.NY), 2],
            domain,
            label = "exact",
        )

        pressure_profile =
            plot(domain, p[x_pos, 1:(problem.NY)], label = "solution", title = "p")
        plot!(
            pressure_profile,
            domain,
            pressure_field[x_pos, 1:(problem.NY)],
            label = "exact",
        )

        temperature_profile =
            plot(domain, T[x_pos, 1:(problem.NY), 1], label = "solution", title = "T")

        sigma_xx_profile = plot(
            domain,
            σ_xx[x_pos, 1:(problem.NY)],
            label = "solution",
            title = "sigma_xx",
        )
        plot!(sigma_xx_profile, domain, σ_xx_field[x_pos, 1:(problem.NY)], label = "exact")

        sigma_xy_profile = plot(
            domain,
            σ_xy[x_pos, 1:(problem.NY)],
            label = "solution",
            title = "sigma_xy",
        )
        plot!(sigma_xy_profile, domain, σ_xy_field[x_pos, 1:(problem.NY)], label = "exact")
    else
        y_pos = round(Int, problem.NY / 2)
        domain = x_range[1:Nx]

        velocity_profile_x = plot(domain, u[:, y_pos, 1], label = "solution", title = "u_x")
        plot!(velocity_profile_x, domain, velocity_field[:, y_pos, 1], label = "exact")

        velocity_profile_y = plot(u[:, y_pos, 2], domain, label = "solution", title = "u_y")
        plot!(velocity_profile_y, velocity_field[:, y_pos, 2], domain, label = "exact")

        pressure_profile = plot(domain, p[:, y_pos], label = "solution", title = "p")
        plot!(pressure_profile, domain, pressure_field[:, y_pos], label = "exact")

        temperature_profile = plot(domain, T[:, y_pos], label = "solution", title = "T")

        sigma_xx_profile =
            plot(domain, σ_xx[:, y_pos], label = "solution", title = "sigma_xx")
        plot!(sigma_xx_profile, domain, σ_xx_field[:, y_pos], label = "exact")

        sigma_xy_profile =
            plot(domain, σ_xy[:, y_pos], label = "solution", title = "sigma_xy")
        plot!(sigma_xy_profile, domain, σ_xy_field[:, y_pos], label = "exact")
    end

    plot(
        contour(ρ, title = "Density", fill = true, cbar = true),
        contour(p, title = "Pressure", fill = true),
        contour(T, title = "Temperature", fill = true),
        # contour(pressure_field', title="pressure analytical", fill=true),
        plot!(streamline(u), title = "Computed"),
        plot!(streamline(velocity_field), title = "exact"),
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
        legend = nothing,
    )
    gui()
end

function streamline(
    j;
    amount_of_arrows = 5,
    step = round(Int, size(j, 1) / amount_of_arrows),
)
    s = (1000, 500)
    velocity_field =
        contour((j[:, :, 1] .^ 2 .+ j[:, :, 2] .^ 2)', fill = true, title = "Momentum")
    N = size(j, 1)
    X = [
        i
        for i in range(2, size(j, 1), step = step), j in range(1, size(j, 2), step = step)
    ]
    Y = [
        j
        for i in range(2, size(j, 1), step = step), j in range(1, size(j, 2), step = step)
    ]

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
