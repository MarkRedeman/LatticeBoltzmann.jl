using LaTeXStrings

struct TakeSnapshots{S <: AbstractVector{Array{Float64, 3}}} <: ProcessingMethod
    problem::FluidFlowProblem
    every_t::Union{Int, Vector{Int}}
    snapshots::S
    timesteps::Vector{Int}
end
TakeSnapshots(problem, every_t) =
    TakeSnapshots(problem, every_t, Array{Float64, 3}[], Int[])

function next!(pm::TakeSnapshots, q, f_in, t::Int64)
    if pm.every_t isa Int && mod(t, pm.every_t) != 0
        return false
    end
    if pm.every_t isa Vector{Int} && t ∉ pm.every_t
        return false
    end

    # snapshot = take_snapshot(pm, q, f_in, t)
    snapshot = copy(f_in)

    push!(pm.snapshots, snapshot)
    push!(pm.timesteps, t)

    visualize(pm, q)

    return false
end

function visualize(pm::TakeSnapshots, q::Quadrature)
    velocity_profile_x = plot(xlabel = L"x", ylabel = L"u_y")
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

    problem = pm.problem
    for (timestep, f_in) in zip(pm.timesteps, pm.snapshots)
        time = delta_t(problem) * timestep
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

            # u[x_idx, y_idx, :] .-= LatticeBoltzmann.velocity(problem, x, y, time)
            σ_ = LatticeBoltzmann.dimensionless_stress(problem, σ_)
            σ_xx[x_idx, y_idx] = σ_[1, 1]
            σ_xy[x_idx, y_idx] = σ_[1, 2]
        end

        x_pos = max(round(Int, problem.NX / 2), 1)
        y_pos = max(round(Int, problem.NY / 2), 1)
        domain = x_range[1:Nx]

        x_range, y_range = range(problem)
        velocity = (x, y, t) -> LatticeBoltzmann.velocity(problem, x, y, t)
        σ = (x, y, t) -> LatticeBoltzmann.deviatoric_tensor(q, problem, x, y, t)
        pr = (x, y, t) -> LatticeBoltzmann.pressure(q, problem, x, y, t)

        exact_range = range(0.0, length = 1000, stop = problem.domain_size[1])

        plot!(
            velocity_profile_x,
            domain,
            u[:, y_pos, 2],
            label = latexstring("t = ", time),
            linecolor = :gray,
        )
        y_exact = map(exact_range) do y
            velocity(y, y_range[y_pos], time)[2]
        end
        # plot!(velocity_profile_x, exact_range, y_exact,  label = "", linecolor = :gray, linealpha = 0.2, linestyle = :dash)

        # annotation_x = 0.5
        # annotation_y = u[x_pos, end - 3, 1]
        # annotation_y = u[x_pos, round(Int, problem.NY /2) + 1, 1] + 0.005
        # annotate!(velocity_profile_x, annotation_y, annotation_x, text(L"\nu / t = 0.5", rotation = 00))

        # plot!(velocity_profile_x, domain, u[x_pos, :, 1], label = latexstring("t = ", time), linecolor = :gray)
        # plot!(velocity_profile_x, (y) -> velocity(x_range[x_pos], y, time)[1], exact_range, label = "", linecolor = :gray, linealpha = 0.2, linestyle = :dash)

        # annotation_x = 0.5
        # annotation_y = u[x_pos, end - 3, 1]
        # annotation_y = u[x_pos, round(Int, problem.NY /2) + 1, 1] + 0.005
        # annotate!(velocity_profile_x, annotation_x, annotation_y, text(L"\nu / t = 0.5", rotation = 00))

        # annotate!(velocity_profile_x, [(5, y[5], Plots.text("this is #5", 16, :red, :center)), (10, y[10], Plots.text("this is #10", :right, 20, "courier"))])

        # scatter!(velocity_profile_x, domain, u[x_pos, :, 2], label = latexstring("t = ", timestep))

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

    ps = (
        velocity_profile_x = velocity_profile_x,
        velocity_profile_y = velocity_profile_y,
        temperature_profile = temperature_profile,
        pressure_profile = pressure_profile,
        sigma_xx_profile = sigma_xx_profile,
        sigma_xy_profile = sigma_xy_profile,
    )

    return ps
end
