module Example
module TaylorGreenVortex

using StatsPlots
using DataFrames
using lbm
using Plots

include("plot.jl")

function process!(tgv, quadrature, f_in, t, stats; visualize = false)
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

    # @inbounds for x in 1:N, y in 1:N
    #     density_field[x, y] = lbm.density(quadrature, tgv, x, y, t)
    #     velocity_field[x, y, :] = lbm.velocity(tgv, x, y, t)
    #     pressure_field[x, y] = lbm.pressure(quadrature, tgv, x, y, t)
    # end

    x_range = range(0, 2pi, length=Nx + 1)
    y_range = range(0, 2pi, length=Ny + 1)
    @inbounds for x_idx in 1:Nx, y_idx in 1:Ny
        x = x_range[x_idx]
        y = y_range[y_idx]

        density_field[x_idx, y_idx] = lbm.density(quadrature, tgv, x, y, t)
        pressure_field[x_idx, y_idx] = lbm.pressure(quadrature, tgv, x, y, t)
        velocity_field[x_idx, y_idx, :] = lbm.velocity(tgv, x, y, t)
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

    # Compare with analytical results?
    push!(stats, [
        sum(ρ),
        sum(j),
        sum(E),
        sum(E_k),
        sum(ϵ),
        sum(density_field),
        sum(velocity_field ./ density_field),
        1.0,
        sum(
            density_field .* (velocity_field[:, :, 1].^2 + velocity_field[:, :, 2].^2)
        ),
        1.0,
        rho_error_squared,
        ux_error_squared,
        uy_error_squared,
        u_error,
    ])

    if visualize == true
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
end

function siumlate(tgv::TaylorGreenVortexExample, quadrature::Quadrature = D2Q9();)
    # NOTE we should base τ on ν and the speed of sound for the given qquadrature

    # initialize
    f_out, collision_operator = lbm.initialize!(quadrature, tgv)
    @show tgv
    @show collision_operator
    f_in = copy(f_out)

    stats = DataFrame(
        [
            Float64[], Float64[], Float64[], Float64[], Float64[],
            Float64[], Float64[], Float64[], Float64[], Float64[],
            Float64[], Float64[], Float64[], Float64[]
        ],
        [
            :density, :momentum, :total_energy, :kinetic_energy, :internal_energy,
            :density_a, :momentum_a, :total_energy_a, :kinetic_energy_a, :internal_energy_a,
            :density_e, :momentum_e, :kinetic_energy_e, :u_error
        ]
    )

    # return f_in
    n_steps = 100 * tgv.scale * tgv.scale
    # n_steps = 10 * tgv.scale * tgv.scale

    @show tgv.ν * (tgv.k_x^2 + tgv.k_y^2)
    @inbounds for t = 0:n_steps
        if mod(t, round(Int, n_steps / 10)) == 0
            @show t, t / n_steps
        end

        if (mod(t, 1) == 0)
            process!(
                tgv,
                quadrature,
                f_in,
                t * tgv.ν * (tgv.k_x^2 + tgv.k_y^2),
                stats,
                visualize = (mod(t, round(Int, n_steps / 100)) == 0)
                # visualize = true
            )
        end

        f_out = collide(collision_operator, quadrature, f_in)

        f_in = stream(quadrature, f_out)

        # check_stability(f_in) || return :unstable, f_in, stats
    end
    process!(tgv, quadrature, f_in, n_steps * tgv.ν * (tgv.k_x^2 + tgv.k_y^2), stats, visualize = true)

    @show stats

    f_in, stats
end

results = let
    stats = DataFrame([Float64[], Int[], Any[]], [:nu, :scale, :stats])

    quadratures = [
        # D2Q4(),
        # D2Q5(),
        D2Q9(),
        # D2Q17(),
    ]

    quadrature = last(quadratures)


    # νs = (0.0:0.5:6.0) ./ 6.0
    # scales = [1, 2, 4, 8]
    # for ν in νs
    # for scale = scales
    # scale = 1
    # ν = 1.0 / 6.0
    # ν = 0.0
    scale = 2
    ν = 1.0 / 6.0
    example = TaylorGreenVortexExample(
        ν,
        scale,
    )

    @time result = siumlate(example, quadrature);
    return result
    push!(stats, [ν, scale, result[2]])
    # end
    # end

    s = stats
    using Plots, LaTeXStrings
    nu_scale_error = Array{Float64}([s.nu s.scale map(i -> i.u_error[end], s.stats)])
    plot()
    for nu_idx = 1:length(νs)
        nu = nu_scale_error[nu_idx * length(scales), 1]
        nu_round = round(nu, digits = 2)

        plot!(
            map(
                i -> 16 * i,
                nu_scale_error[(nu_idx - 1) * length(scales) .+ (1:length(scales)), 2]
            ),
            nu_scale_error[(nu_idx - 1) * length(scales) .+ (1:length(scales)), 3],
            label="nu = $nu_round"
        )
    end
    plot!(x -> 10.0 * x^(-2), label="y(x) = 10x^(-2)", linestyle=:dash)
    plot!(
        scale=:log10,
        legend=:bottomleft,
        legendfontsize=5
    )
    gui()

    nu_idx = 0
    plot(nu_scale_error[nu_idx * length(scales) .+ (1:length(scales)), 2], nu_scale_error[nu_idx * length(scales) .+ (1:length(scales)), 3])
end

# Some thoughs: hide the storage of the distributions f inside of an interface
# so that we can do:
# for x in xs, y in ys, z in zs
#  f = fs[x, y ,z]::Vector
#
# the interface could hide it by returning
# return @view f_internal[x, y, z, :]


# Idea: introduce an Initial Value Problem
# problem = InitialValueProblem
# solution = solve(problem, LBM(Lattice, CollisionModel))
# LBM(Lattice, CollisionModel) can be a solution method

end
end
