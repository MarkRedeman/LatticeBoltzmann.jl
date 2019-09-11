module Example
module TaylorGreenVortex

using StatsPlots
using DataFrames
using lbm
using Plots

struct TaylorGreenVortexExample <: lbm.InitialValueProblem
    scale
    rho_0
    u_max
    ν
    NX
    NY
    k_x
    k_y

    function TaylorGreenVortexExample(ν = 1.0 / 6.0 , scale = 2, NX = 16 * scale, NY = NX)
        u_max = 0.02 / scale
        Re = NX * u_max / ν
        @show Re
        return new(
            scale,
            1.0,
            u_max,
            ν,
            NX,
            NY,
            2pi / NX,
            2pi / NY,
        )
    end
end
function density(tgv::TaylorGreenVortexExample, x::Int, y::Int, timestep::Int = 0)
    X = x + 0.5
    Y = y + 0.5
    u_max = tgv.u_max
    kx = tgv.k_x
    ky = tgv.k_y

    P = -0.25 * tgv.rho_0 * u_max * u_max * (
        (ky / kx) * cos(2.0 * kx * X) + (kx / ky)*cos(2.0 * ky * Y)
    ) * decay(tgv, x, y, timestep)^2;

    return 1.0 + 3.0 * P
end
function velocity(tgv::TaylorGreenVortexExample, x::Int, y::Int, timestep::Int = 0)
    X = x + 0.5
    Y = y + 0.5
    u_max = tgv.u_max
    kx = tgv.k_x
    ky = tgv.k_y

    return decay(tgv, x, y, timestep) .* [
      -u_max * sqrt(ky / kx) * cos(kx * X) * sin(ky * Y),
       u_max * sqrt(kx / ky) * sin(kx * X) * cos(ky * Y)
    ]
end
function decay(tgv::TaylorGreenVortexExample, x, y, timestep::Int)
    td = 1.0 / (tgv.ν * (tgv.k_x^2 + tgv.k_y^2));

    exp(-1.0 * timestep / td)
end


function process!(tgv, quadrature, f_in, t, stats)
    # Density
    ρ = lbm.density(quadrature, f_in)

    # Momentum
    j = lbm.momentum(quadrature, f_in)
    E = lbm.total_energy(quadrature, f_in)
    E_k = lbm.kinetic_energy(quadrature, f_in, ρ, j ./ ρ)
    ϵ = 1.0 #lbm.internal_energy(quadrature, f_in, ρ, j ./ ρ)

    # E_k = sum((j[:, :, 1].^2 + j[:, :, 1].^2) ./ ρ)


    N = size(f_in, 1)
    density_field = fill(0.0, N, N)
    velocity_field = fill(0.0, N, N, lbm.dimension(quadrature))

    for x in 1:N, y in 1:N
        density_field[x, y] = density(tgv, x, y, t)
        velocity_field[x, y, :] = velocity(tgv, x, y, t)
    end

    rho_error_squared = sqrt(
        sum(
            (ρ - density_field).^2 ./
            (density_field).^2
        )
    )
    ux_error_squared = sqrt(
        sum(
            (ρ .* j[:, :, 1] - velocity_field[:, :, 1]).^2 ./ velocity_field[:, :, 1].^2
        )
    )
    uy_error_squared = sqrt(
        sum(
            (ρ .* j[:, :, 2] - velocity_field[:, :, 2]).^2 ./ velocity_field[:, :, 2].^2
        )
    )
    u_error = sqrt(
        sum(
            ((ρ .* j[:, :, 2] - velocity_field[:, :, 2]).^2 .+
             (ρ .* j[:, :, 2] - velocity_field[:, :, 2]).^2) #./
            # (velocity_field[:, :, 1].^2 .+ velocity_field[:, :, 2].^2)
        )
    )


    # sumrhoe2 += (rho - rhoa) * (rho - rhoa);
    # sumuxe2 += (ux - uxa) * (ux - uxa);
    # sumuye2 += (uy - uya) * (uy - uya);

    # sumrhoa2 += (rhoa - rho0) * (rhoa - rho0);
    # sumuxa2 += uxa * uxa;
    # sumuya2 += uya * uya;

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

    # if (mod(t, 50) == 0)
    #     plot(
    #         plot(stats.density),
    #         plot(stats.momentum),
    #         # plot(stats.total_energy),
    #         plot(stats.kinetic_energy),
    #         plot(stats.density_a),
    #         plot(stats.momentum_a),
    #         # plot(stats.total_energy),
    #         plot(stats.kinetic_energy_a),
    #         # plot(stats.internal_energy),
    #         legend=nothing
    #     )

    #     gui()
    # end

    u = j[:, :, 1] ./ ρ
    v = j[:, :, 2] ./ ρ

    s = (1000, 500)

    domain = (1 : size(j, 1)) ./ size(j, 1)
    velocity_profile_x = plot(domain, j[:, 4, 1], size=s)
    plot!(velocity_profile_x, domain, velocity_field[:, 4, 1], size=s)
    velocity_profile_y = plot(j[:, 4, 2], domain, size=s)
    plot!(velocity_profile_y, velocity_field[:, 4, 2], domain, size=s)


    plot(
        contour(ρ[:, :, 1], fill=true, clims=(0, 1.05), cbar=true, size=s),
        streamline(j),
        streamline(velocity_field),
        streamline(velocity_field .- ρ .* j),
        # velocity_field,
        plot(stats.u_error, legend=false, title="U_e"),
        plot(stats.kinetic_energy, legend=false, title="Kinetic energy", size=s),
        # plot(stats.total_energy, legend=false, title="total energy", size=s),
        # contour(u[:, :, 1].^2, fill=true, cbar=true, size=s, title="u"),
        # contour(v[:, :, 1].^2, fill=true, cbar=true, size=s, title="v"),
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
            # quiver=(u[:, :, 1], v[:, :, 1])
            quiver=(x, y) -> (j[x, y, 1] / sqrt(j[x, y, 1]^2 + j[x, y, 2]^2), j[x, y, 2] / sqrt(j[x, y, 1]^2 + j[x, y, 2]^2)) ,
            color="white",
        )
    return velocity_field
end

function initialize!(quadrature, tgv, N)
    density_field = fill(0.0, N, N)
    velocity_field = fill(0.0, N, N, lbm.dimension(quadrature))
    for x in 1:N, y in 1:N
        density_field[x, y] = density(tgv, x, y)
        velocity_field[x, y, :] = velocity(tgv, x, y)
    end

    return lbm.equilibrium(
        quadrature,
        density_field,
        velocity_field,
        1.0
    );
    # Add offequilibrium ?
end

function siumlate(tgv::TaylorGreenVortexExample, quadrature::Quadrature = D2Q9();)
    # NOTE we should base τ on ν and the speed of sound for the given qquadrature
    @show tgv
    N = tgv.NX
    τ = quadrature.speed_of_sound_squared * tgv.ν + 0.5
    @show τ
    # τ = 0.5;

    # initialize
    f_out = initialize!(quadrature, tgv, N)
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
    @inbounds for t = 0 : n_steps
        if mod(t, round(Int, n_steps / 10)) == 0
            @show t, t / n_steps
        end

        if (mod(t, 1) == 10)
            process!(tgv, quadrature, f_in, t, stats)
        end

        f_out = collide(SRT(τ), quadrature, f_in)

        f_in = stream(quadrature, f_out)
    end
    process!(tgv, quadrature, f_in, n_steps + 1, stats)

    @show stats

    f_in, stats
end

stats = DataFrame([Float64[], Int[], Any[]], [:nu, :scale, :stats])

quadratures = [
    D2Q4(),
    D2Q5(),
    D2Q9(),
    D2Q17(),
]

quadrature = last(quadratures)


νs = (0.0:6.0) ./ 6.0
scales = [1, 2, 4]
for ν in νs
for scale = scales
# scale = 1
# ν = 1.0 / 6.0
# ν = 0.0
    example = TaylorGreenVortexExample(
        ν,
        scale,
    )

    @time result = siumlate(example, quadrature);
    push!(stats, [ν, scale, result[2]])
end
end

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
# plot!(
#     map(
#         i -> 16 * i,
#         scales
#     ),
#     map(i -> 1E-1 * 4.0^(-i), 1:length(scales)),
#     label="4^(-x)"
# )
gui()

nu_idx = 0
    plot(nu_scale_error[nu_idx * length(scales) .+ (1:length(scales)), 2], nu_scale_error[nu_idx * length(scales) .+ (1:length(scales)), 3])

end
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
