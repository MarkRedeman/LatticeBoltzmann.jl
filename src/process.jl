using Plots
using DataFrames

abstract type ProcessingMethod end
struct CompareWithAnalyticalSolution5{T} <: ProcessingMethod
    problem::InitialValueProblem
    should_process::Bool
    n_steps::Int64
    stop_criteria::StopCriteria
    df::T
end
CompareWithAnalyticalSolution(problem, should_process, n_steps) = CompareWithAnalyticalSolution5(
    problem,
    should_process,
    n_steps,
    StopCriteria(problem),
    DataFrame(
        [
            Float64[], Float64[], Float64[], Float64[], Float64[],
            Float64[], Float64[], Float64[], Float64[], Float64[],
            Float64[], Float64[]
        ],
        [
            :density, :momentum, :total_energy, :kinetic_energy, :internal_energy,
            :density_a, :momentum_a, :total_energy_a, :kinetic_energy_a, :internal_energy_a,
            :u_error, :p_error
        ]
    )
)

function next!(stats::CompareWithAnalyticalSolution5, q, f_in, t::Int64)
    if mod(t, 100) == 0
        if (should_stop!(stats.stop_criteria, q, f_in))
            @info "Stopping after $t steps out of $n_steps"
            return true
        end
    end

    if (! stats.should_process)
        return false
    end

    Δt = delta_t(stats.problem)
    time = t * Δt

    if mod(t, 1) == 0
        process!(
            stats.problem,
            q,
            f_in,
            time,
            stats.df;
            should_visualize = (
                t == stats.n_steps ||
                mod(t, round(Int, stats.n_steps / 20)) == 0
            )
        )
    end

    return false
end

function process!(problem::InitialValueProblem, q::Quadrature, f_in, time, stats; should_visualize = false)
    f = Array{Float64}(undef, size(f_in, 3))
    u = zeros(dimension(q))
    expected_u = zeros(dimension(q))

    Nx = size(f_in, 1)
    Ny = size(f_in, 2)

    x_range, y_range = range(problem)

    total_density = 0.0
    total_momentum = 0.0
    total_energy = 0.0
    total_kinetic_energy = 0.0
    total_internal_energy = 0.0

    expected_total_density = 0.0
    expected_total_momentum = 0.0
    expected_total_energy = 0.0
    expected_total_kinetic_energy = 0.0
    expected_total_internal_energy = 0.0

    rho_error_squared = 0.0
    ux_error_squared = 0.0
    uy_error_squared = 0.0
    u_error = 0.0
    p_error = 0.0

    @inbounds for x_idx in 1:Nx, y_idx in 1:Ny
        # Calculated
        @inbounds for f_idx = 1 : size(f_in, 3)
            f[f_idx] = f_in[x_idx, y_idx, f_idx]
        end
        ρ = density(q, f)
        velocity!(q, f, ρ, u)

        # Adding the forcing term moves the optimal tau for poiseuille flows
        # F = cm.force(x_idx, y_idx, 0.0)
        # u += cm.τ * F

        T = temperature(q, f, ρ, u)
        p = pressure(q, f, ρ, u)

        ρ = dimensionless_density(problem, ρ)
        u = dimensionless_velocity(problem, u)
        T = dimensionless_temperature(q, problem, T)
        p = dimensionless_pressure(q, problem, p)

        total_density += ρ
        total_momentum += (u[1] + u[2]) * ρ
        kinetic_energy = (u[1]^2 + u[2]^2) * ρ
        internal_energy = T

        total_kinetic_energy += kinetic_energy
        total_internal_energy += internal_energy
        total_energy += kinetic_energy + internal_energy

        # Analytical
        x = x_range[x_idx]
        y = y_range[y_idx]

        # Compute statistics
        expected_ρ = lbm.density(q, problem, x, y, time)
        expected_p = lbm.pressure(q, problem, x, y, time)
        expected_ϵ = (dimension(q) / 2) * expected_p / expected_ρ
        expected_v = lbm.velocity(problem, x, y, time)
        expected_T = expected_p / expected_ρ

        # @show u[1], expected_v[1]

        expected_kinetic_energy = (expected_v[1]^2 + expected_v[2]^2)
        expected_internal_energy = expected_T

        expected_total_density += expected_ρ
        expected_total_momentum += expected_ρ * (expected_v[1] + expected_v[2])
        expected_total_kinetic_energy += expected_kinetic_energy
        expected_total_internal_energy += expected_internal_energy
        expected_total_energy += expected_kinetic_energy + expected_internal_energy

        rho_error_squared += (ρ - expected_ρ)^2
        ux_error_squared += (u[1] - expected_v[1])^2
        uy_error_squared += (u[2] - expected_v[2])^2
        # u_error += sqrt((u[1] - expected_v[1])^2 + (u[2] - expected_v[2])^2)
        opp = Float64(y_range.step) * Float64(x_range.step)
        # opp = 1.0
        # opp = Float64(y_range.step)
        u_error += (
            opp * (
                (
                    (u[1] - expected_v[1])
                )^2 +
                (
                    (u[2] - expected_v[2])
                )^2
            )
        )
        p_error += opp * (
            (p - expected_p)^2
        )
    end

    # Compare with analytical results?
    push!(stats, [
        total_density,
        total_momentum,
        total_energy,
        total_kinetic_energy,
        total_internal_energy,
        expected_total_density,
        expected_total_momentum,
        expected_total_energy,
        expected_total_kinetic_energy,
        expected_total_internal_energy,
        # rho_error_squared,
        # ux_error_squared,
        # uy_error_squared,
        # u_error / expected_total_momentum,
        # sqrt(u_error) / expected_total_momentum,
        sqrt(u_error), #/ expected_total_momentum
        sqrt(p_error)
        # sqrt(ux_error_squared),
        # sqrt(uy_error_squared)
    ])

    if should_visualize
        visualize(problem, q, f_in, time, stats)
    end
    return
end

function visualize(problem::InitialValueProblem, quadrature::Quadrature, f_in, time, stats)
    q = quadrature
    # Density
    ρ = lbm.density(quadrature, f_in)
    ρ = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))
    p = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))
    T = Array{Float64}(undef, size(f_in, 1), size(f_in, 2))

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

    x_range, y_range = range(problem)

    f = Array{Float64}(undef, size(f_in, 3))
    u = zeros(dimension(q))
    @inbounds for x_idx in 1:Nx, y_idx in 1:Ny
        @inbounds for f_idx = 1 : size(f_in, 3)
            f[f_idx] = f_in[x_idx, y_idx, f_idx]
        end
        x = x_range[x_idx]
        y = y_range[y_idx]

        density_field[x_idx, y_idx] = lbm.density(quadrature, problem, x, y, time)
        pressure_field[x_idx, y_idx] = lbm.pressure(quadrature, problem, x, y, time)
        velocity_field[x_idx, y_idx, :] = lbm.velocity(problem, x, y, time)

        ρ[x_idx, y_idx] = density(q, f)
        p[x_idx, y_idx] = pressure(q, f, ρ[x_idx, y_idx], j[x_idx, y_idx, :] ./ ρ[x_idx, y_idx])
        T[x_idx, y_idx] = temperature(q, f, ρ[x_idx, y_idx], j[x_idx, y_idx, :] ./ ρ[x_idx, y_idx])
        # velocity!(q, f, ρ, u)

        ρ[x_idx, y_idx] = dimensionless_density(problem, ρ[x_idx, y_idx])
        p[x_idx, y_idx] = dimensionless_pressure(q, problem, p[x_idx, y_idx])
        T[x_idx, y_idx] = dimensionless_temperature(q, problem, T[x_idx, y_idx])
        # T[x_idx, y_idx] = dimensionless_temperature(problem, T[x_idx, y_idx])
        j[x_idx, y_idx, :] = dimensionless_velocity(problem, j[x_idx, y_idx, :] ./ ρ[x_idx, y_idx])
    end

    s = (1000, 500)

    domain = (2 : (problem.NY - 1)) ./ (problem.NY - 2)

    if (typeof(problem) != DecayingShearFlow)
        x_pos = round(Int, problem.NX / 2)
        domain = y_range[1:Ny]

        velocity_profile_x = plot(domain, j[x_pos, 1:(problem.NY), 1], size=s, label="solution", title="u_x", legend=nothing)
        plot!(velocity_profile_x, domain, velocity_field[x_pos, 1:(problem.NY), 1], size=s, label="exact")

        velocity_profile_y = plot(j[x_pos, 1:(problem.NY), 2], domain, size=s, label="solution", title="u_y", legend=nothing)
        plot!(velocity_profile_y, velocity_field[x_pos, 1:(problem.NY), 2], domain, size=s, label="exact")

        pressure_profile = plot(domain, p[x_pos, 1:(problem.NY)], size=s, label="solution", title="p", legend=nothing)
        plot!(pressure_profile, domain, pressure_field[x_pos, 1:(problem.NY)], size=s, label="exact")

        temperature_profile = plot(domain, T[x_pos, 1:(problem.NY), 1], size=s, label="solution", title="T", legend=nothing)
        # plot!(temperature_profile, domain, temperature_field[x_pos, 1:(problem.NY), 1], size=s, label="exact")
    else
        y_pos = round(Int, problem.NY / 2)
        domain = x_range[1:Nx]

        velocity_profile_x = plot(domain, j[:, y_pos, 1], size=s, label="solution", title="u_x")
        plot!(velocity_profile_x, domain, velocity_field[:, y_pos, 1], size=s, label="exact")

        velocity_profile_y = plot(j[:, y_pos, 2], domain, size=s, label="solution", title="u_y")
        plot!(velocity_profile_y, velocity_field[:, y_pos, 2], domain, size=s, label="exact")

        pressure_profile = plot(domain, p[:, y_pos], size=s, label="solution", title="p", legend=nothing)
        plot!(pressure_profile, domain, pressure_field[:, y_pos], size=s, label="exact")

        temperature_profile = plot(domain, T[:, y_pos], size=s, label="solution", title="T", legend=nothing)
    end


    kinetic_energy_profile = plot(stats.kinetic_energy, legend=false, title="Kinetic energy", size=s)

    plot(
        contour(ρ, title="Density", fill=true, cbar=true, size=s),
        contour(p, title="Pressure", fill=true),
        contour(T, title="Temperature", fill=true),
        # contour(pressure_field', title="pressure analytical", fill=true),
        plot!(streamline(j), title="Computed"),
        plot!(streamline(velocity_field), title="exact"),
        # heatmap(j[:, :, 1]', fill=true),
        # heatmap(j[:, 2:(problem.NY-1), 2]', fill=true),
        # heatmap(velocity_field[:, :, 1]', fill=true),
        # heatmap(velocity_field[:, :, 2]', fill=true),
        # streamline(velocity_field .- ρ .* j),
        plot(stats.u_error, legend=false, title="U_e"),
        plot(stats.p_error, legend=false, title="P_e"),
        kinetic_energy_profile,
        pressure_profile,
        temperature_profile,
        velocity_profile_x,
        velocity_profile_y,
        size=(1000, 600)
    )
    gui()
end

function streamline(j; amount_of_arrows = 5, step = round(Int, size(j, 1) / amount_of_arrows) )
    s = (1000, 500)
    velocity_field = contour((j[:, :, 1].^2 .+ j[:, :, 2].^2)', cbar=true, fill=true, title="Momentum", size=s)
    N = size(j, 1)
    X = [i for i in range(2, size(j, 1), step = step), j in range(1, size(j, 2), step = step)]
    Y = [j for i in range(2, size(j, 1), step = step), j in range(1, size(j, 2), step = step)]
    # @show "process: ", u, v,

    quiver!(
        velocity_field,
        X, Y,
        quiver=(x, y) -> (
            0.1 * amount_of_arrows * j[Int(x), Int(y), 1] / sqrt(j[Int(x), Int(y), 1]^2 + j[Int(x), Int(y), 2]^2),
            0.1 * amount_of_arrows * j[Int(x), Int(y), 2] / sqrt(j[Int(x), Int(y), 1]^2 + j[Int(x), Int(y), 2]^2)
        ),
        color="white",
    )

    return velocity_field
end
