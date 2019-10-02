using Plots

export process!, initialize, apply_boundary_conditions!,
    InitialValueProblem, viscosity, delta_t,
    lattice_velocity,
    lattice_density,
    lattice_pressure,
    lattice_force,
    lattice_viscosity

abstract type InitialValueProblem end

has_external_force(problem::InitialValueProblem) = false

function initialize(quadrature::Quadrature, lattice::Lattice, problem::InitialValueProblem)

end

function initial_equilibrium(quadrature::Quadrature, problem::InitialValueProblem, x::Float64, y::Float64)
    return equilibrium(
        quadrature,
        lattice_density(quadrature, problem, x, y),
        lattice_velocity(quadrature, problem, x, y),
        lattice_pressure(quadrature, problem, x, y) / lattice_density(quadrature, problem, x, y)
    )
end

function initial_condition(q::Quadrature, problem::InitialValueProblem, x::Float64, y::Float64)
    initial_equilibrium(q, problem, x, y)
end

function initialize(quadrature::Quadrature, problem::InitialValueProblem)
    # force_field = Array{Any}(undef, problem.NX, problem.NY)
    f = Array{Float64}(undef, problem.NX, problem.NY, length(quadrature.weights))

    # NOTE: we have periodic boundaries
    x_range = range(0, problem.domain_size[1], length=problem.NX + 1)
    y_range = range(0, problem.domain_size[2], length=problem.NY + 1)
    for x_idx in 1:problem.NX, y_idx in 1:problem.NY
        f[x_idx, y_idx, :] = initial_equilibrium(
            quadrature,
            problem,
            x_range[x_idx],
            y_range[y_idx]
        )

        # force_field[x_idx, y_idx, :] = force(problem, x, y)
    end
    force_field = (x_idx, y_idx, t) -> lattice_force(problem, x_idx, y_idx, t)

    τ = quadrature.speed_of_sound_squared * lattice_viscosity(problem) + 0.5
    @show τ
    if has_external_force(problem)
        collision_operator = SRT_Force(τ, force_field)
    else
        collision_operator = SRT(τ)
    end

    return f, collision_operator
end

function delta_t(problem)
    ν = viscosity(problem)
    k_x = problem.domain_size[1] / problem.NX
    k_y = problem.domain_size[2] / problem.NY

    Δt = ν * (k_x^2 + k_y^2)

    return Δt
end

function process!(problem::InitialValueProblem, q::Quadrature, f_in, time, stats; should_visualize = false)
    f = Array{Float64}(undef, size(f_in, 3))
    u = zeros(dimension(q))
    expected_u = zeros(dimension(q))

    Nx = size(f_in, 1)
    Ny = size(f_in, 2)

    x_range = range(0, problem.domain_size[1], length=Nx + 1)
    y_range = range(0, problem.domain_size[2], length=Ny + 1)

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

    @inbounds for x_idx in 1:Nx, y_idx in 1:Ny
        # Calculated
        @inbounds for f_idx = 1 : size(f_in, 3)
            f[f_idx] = f_in[x_idx, y_idx, f_idx]
        end
        ρ = density(q, f)
        velocity!(q, f, ρ, u)
        T = temperature(q, f, ρ, u)

        ρ = dimensionless_density(problem, ρ)
        u = dimensionless_velocity(problem, u)
        T = dimensionless_temperature(problem, T)

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
                    (u[1] - expected_v[1]) / problem.u_max
                )^2 +
                (
                    (u[2] - expected_v[2]) / problem.u_max
                )^2
            )
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
        sqrt(u_error) #/ expected_total_momentum
        # sqrt(ux_error_squared),
        # sqrt(uy_error_squared)
    ])

    if should_visualize
        visualize(problem, q, f_in, time, stats)
    end
    return
end

function visualize(problem::InitialValueProblem, quadrature::Quadrature, f_in, time, stats)
    # Density
    ρ = lbm.density(quadrature, f_in)

    # Momentum
    j = lbm.momentum(quadrature, f_in)
    E = lbm.total_energy(quadrature, f_in)
    E_k = lbm.kinetic_energy(quadrature, f_in, ρ, j ./ ρ)
    ϵ = 1.0 #lbm.internal_energy(quadrature, f_in, ρ, j ./ ρ)


    Nx = size(f_in, 1)
    Ny = size(f_in, 2)
    density_field = fill(0.0, Nx, Ny)
    pressure_field = fill(0.0, Nx, Ny)
    velocity_field = fill(0.0, Nx, Ny, lbm.dimension(quadrature))

    x_range = range(0, problem.domain_size[1], length=Nx + 1)
    y_range = range(0, problem.domain_size[2], length=Ny + 1)
    @inbounds for x_idx in 1:Nx, y_idx in 1:Ny
        x = x_range[x_idx]
        y = y_range[y_idx]

        density_field[x_idx, y_idx] = lbm.density(quadrature, problem, x, y, time)
        pressure_field[x_idx, y_idx] = lbm.pressure(quadrature, problem, x, y, time)
        velocity_field[x_idx, y_idx, :] = lbm.velocity(problem, x, y, time)

        ρ[x_idx, y_idx] = dimensionless_density(problem, ρ[x_idx, y_idx])
        j[x_idx, y_idx, :] = dimensionless_velocity(problem, j[x_idx, y_idx, :] ./ ρ[x_idx, y_idx])
    end

    s = (1000, 500)

    domain = (1 : problem.NX) ./ problem.NX
    velocity_profile_x = plot(domain, j[:, 4, 1], size=s, label="solution")
    plot!(velocity_profile_x, domain, velocity_field[:, 4, 1], size=s, label="exact")
    velocity_profile_y = plot(j[:, 4, 2], domain, size=s, label="solution")
    plot!(velocity_profile_y, velocity_field[:, 4, 2], domain, size=s, label="exact")

    kinetic_energy_profile = plot(stats.kinetic_energy, legend=false, title="Kinetic energy", size=s)

    plot(
        # contour(ρ[:, :, 1], fill=true, clims=(0, 1.05), cbar=true, size=s),
        # contour(lbm.pressure(quadrature, f_in, ρ[:, :, 1], j), title="pressure"),
        # contour(pressure_field, title="pressure analytical"),
        plot!(streamline(j), title="Solution"),
        plot!(streamline(velocity_field), title="exact"),
        # streamline(velocity_field .- ρ .* j),
        plot(stats.u_error, legend=false, title="U_e"),
        kinetic_energy_profile,
        velocity_profile_x,
        velocity_profile_y,
        size=(1000, 600)
    )
    gui()
end

function streamline(j; amount_of_arrows = 10, step = round(Int, size(j, 1) / amount_of_arrows) )
    s = (1000, 500)
    velocity_field = contour(j[:, :, 1].^2 .+ j[:, :, 2].^2, cbar=true, fill=true, title="Momentum", size=s)
    N = size(j, 1)
    X = [i for i in range(1, size(j, 1), step = step), j in range(1, size(j, 2), step = step)]
    Y = [j for i in range(1, size(j, 1), step = step), j in range(1, size(j, 2), step = step)]
    # @show "process: ", u, v,

    quiver!(
        velocity_field,
        X, Y,
        quiver=(x, y) -> (
            0.5 * amount_of_arrows * j[x, y, 1] / sqrt(j[x, y, 1]^2 + j[x, y, 2]^2),
            0.5 * amount_of_arrows * j[x, y, 2] / sqrt(j[x, y, 1]^2 + j[x, y, 2]^2)
        ),
        color="white",
    )

    return velocity_field
end

function is_fluid(problem::InitialValueProblem, x, y)
    return true
end

function apply_boundary_conditions!(q::Quadrature, problem::InitialValueProblem, f_in, f_out; time = 0.0)
    nothing
end

function force(problem::InitialValueProblem, x_idx::Int64, y_idx::Int64, time::Float64 = 0.0)
    x_range = range(0, problem.domain_size[1], length=problem.NX + 1)
    y_range = range(0, problem.domain_size[2], length=problem.NY + 1)

    x = x_range[x_idx]
    y = y_range[y_idx]

    return force(problem, x, y, time)
end

function force(problem::InitialValueProblem, x_idx::Int64, y_idx::Int64, time::Int64 = 0.0)
    x_range = range(0, problem.domain_size[1], length=problem.NX + 1)
    y_range = range(0, problem.domain_size[2], length=problem.NY + 1)
    Δt = lbm.delta_t(problem)

    x = x_range[x_idx]
    y = y_range[y_idx]

    return force(problem, x, y, Δt * time)
end

# Dimensionless
function viscosity(problem) #::InitialValueProblem)
    return problem.ν * delta_x(problem)^2 / delta_t(problem)
end

function delta_t(problem::InitialValueProblem)
    return delta_x(problem) * problem.u_max
end

function delta_x(problem::InitialValueProblem)
    return problem.domain_size[1] * (1 / problem.NX)
end

lattice_viscosity(problem) = problem.ν #::InitialValueProblem)
lattice_density(q, problem::InitialValueProblem, x, y, t = 0.0) = density(q, problem, x, y, t)
lattice_velocity(q, problem::InitialValueProblem, x, y, t = 0.0) = problem.u_max * velocity(problem, x, y, t)
lattice_pressure(q, problem::InitialValueProblem, x, y, t = 0.0) = pressure(q, problem, x, y, t)
lattice_force(problem::InitialValueProblem, x, y, t = 0.0) = problem.u_max * delta_t(problem) * force(problem, x, y, t)

dimensionless_viscosity(problem) = problem.ν * delta_x(problem)^2 / delta_t(problem)
dimensionless_density(problem::InitialValueProblem, ρ) = ρ
dimensionless_velocity(problem::InitialValueProblem, u) = u / problem.u_max
dimensionless_pressure(problem::InitialValueProblem, p) = p
dimensionless_temperature(problem::InitialValueProblem, T) = T
dimensionless_force(problem::InitialValueProblem, F) = F / (problem.u_max * delta_t(problem))

include("taylor-green-vortex-decay.jl")
include("decaying-shear-flow.jl")

# error(::Val{:density}, node, solution) = density(node) - density(solution)
