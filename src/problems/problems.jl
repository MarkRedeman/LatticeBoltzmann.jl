using Plots

export process!, initialize, InitialValueProblem, viscosity, delta_t

abstract type InitialValueProblem end
function initialize(quadrature::Quadrature, lattice::Lattice, problem::InitialValueProblem)

end

function viscosity(problem)
    return problem.ν
end

function delta_t(problem)
    ν = viscosity(problem)
    Δt = ν * (problem.k_x^2 + problem.k_y^2)

    return Δt
end

function process!(tgv::InitialValueProblem, quadrature::Quadrature, f_in, time, stats; should_visualize = false)
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

    x_range = range(0, 2pi, length=Nx + 1)
    y_range = range(0, 2pi, length=Ny + 1)
    @inbounds for x_idx in 1:Nx, y_idx in 1:Ny
        x = x_range[x_idx]
        y = y_range[y_idx]

        density_field[x_idx, y_idx] = lbm.density(quadrature, tgv, x, y, time)
        pressure_field[x_idx, y_idx] = lbm.pressure(quadrature, tgv, x, y, time)
        velocity_field[x_idx, y_idx, :] = lbm.velocity(tgv, x, y, time)
    end

    rho_error_squared = sqrt(
        sum((ρ - density_field).^2 ./ (density_field).^2)
    )
    ux_error_squared = sqrt(
        sum((ρ .* j[:, :, 1] - velocity_field[:, :, 1]).^2 ./ velocity_field[:, :, 1].^2)
    )
    uy_error_squared = sqrt(
        sum((ρ .* j[:, :, 2] - velocity_field[:, :, 2]).^2 ./ velocity_field[:, :, 2].^2)
    )
    u_error = sqrt(
        sum(
            ((ρ .* j[:, :, 1] - velocity_field[:, :, 1]).^2 .+
             (ρ .* j[:, :, 2] - velocity_field[:, :, 2]).^2) #./
            # (velocity_field[:, :, 1].^2 .+ velocity_field[:, :, 2].^2)
        )
    )

    total_density = sum(ρ)
    total_momentum = sum(j)
    total_energy = sum(E)
    total_kinetic_energy = sum(E_k)
    total_internal_energy = sum(ϵ)

    expected_total_density = sum(density_field)
    expected_total_momentum = sum(velocity_field ./ density_field)
    expected_total_energy = 1.0
    expected_total_kinetic_energy = sum(density_field .* (velocity_field[:, :, 1].^2 + velocity_field[:, :, 2].^2))
    expected_total_internal_energy = 1.0

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
        rho_error_squared,
        ux_error_squared,
        uy_error_squared,
        u_error,
    ])


    if should_visualize == true
        visualize(tgv, quadrature, f_in, time, stats)
    end
end

function visualize(tgv::InitialValueProblem, quadrature::Quadrature, f_in, time, stats)
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

    x_range = range(0, 2pi, length=Nx + 1)
    y_range = range(0, 2pi, length=Ny + 1)
    @inbounds for x_idx in 1:Nx, y_idx in 1:Ny
        x = x_range[x_idx]
        y = y_range[y_idx]

        density_field[x_idx, y_idx] = lbm.density(quadrature, tgv, x, y, time)
        pressure_field[x_idx, y_idx] = lbm.pressure(quadrature, tgv, x, y, time)
        velocity_field[x_idx, y_idx, :] = lbm.velocity(tgv, x, y, time)
    end

    s = (1000, 500)

    domain = (1 : size(j, 1)) ./ size(j, 1)
    velocity_profile_x = plot(domain, j[:, 4, 1], size=s)
    plot!(velocity_profile_x, domain, velocity_field[:, 4, 1], size=s)
    velocity_profile_y = plot(j[:, 4, 2], domain, size=s)
    plot!(velocity_profile_y, velocity_field[:, 4, 2], domain, size=s)

    kinetic_energy_profile = plot(stats.kinetic_energy, legend=false, title="Kinetic energy", size=s)

    plot(
        contour(ρ[:, :, 1], fill=true, clims=(0, 1.05), cbar=true, size=s),
        contour(lbm.pressure(quadrature, f_in, ρ[:, :, 1], j), title="pressure"),
        contour(pressure_field, title="pressure analytical"),
        streamline(j),
        streamline(velocity_field),
        streamline(velocity_field .- ρ .* j),
        plot(stats.u_error, legend=false, title="U_e"),
        kinetic_energy_profile,
        velocity_profile_x,
        velocity_profile_y,
        size=(1000, 600)
    )
    gui()
end

function streamline(j; step = round(Int, size(j, 1) / 5) )
    s = (1000, 500)
    velocity_field = contour(j[:, :, 1].^2 .+ j[:, :, 2].^2, cbar=true, fill=true, title="Momentum", size=s)
    N = size(j, 1)
    X = [i for i in range(1, size(j, 1), step = step), j in range(1, size(j, 2), step = step)]
    Y = [j for i in range(1, size(j, 1), step = step), j in range(1, size(j, 2), step = step)]
    # @show "process: ", u, v,

    quiver!(
        velocity_field,
        X, Y,
        quiver=(x, y) -> (j[x, y, 1] / sqrt(j[x, y, 1]^2 + j[x, y, 2]^2), j[x, y, 2] / sqrt(j[x, y, 1]^2 + j[x, y, 2]^2)) ,
        color="white",
    )

    return velocity_field
end

include("taylor-green-vortex-decay.jl")

# error(::Val{:density}, node, solution) = density(node) - density(solution)
