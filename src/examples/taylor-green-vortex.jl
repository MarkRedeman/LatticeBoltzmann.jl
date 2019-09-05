module Example
module TaylorGreenVortex

using StatsPlots
using DataFrames
using lbm
using Plots

struct TaylorGreenVortexExample <: lbm.InitialValueProblem
    rho_0
    u_max
    ν
    NX
    NY
    k_x
    k_y

    function TaylorGreenVortexExample(ν = 1.0 / 6.0 , scale = 2, NX = 32 * scale, NY = NX)
        u_max = 0.04 / scale
        Re = NX * u_max / ν
        @show Re
        return new(
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

function initialize!(quadrature, tgv, N)
    density_field = fill(0.0, N, N)
    velocity_field = fill(0.0, N, N, dimension(quadrature))
    for x in 1:N, y in 1:N
        density_field[x, y] = density(tgv, x, y)
        velocity_field[x, y, :] = velocity(tgv, x, y)
    end

    return equilibrium(
        quadrature,
        density_field,
        velocity_field,
        1.0
    );
    # Add offequilibrium ?
end

function siumlate(tgv::TaylorGreenVortexExample;)
    quadratures = [
        D2Q4(),
        D2Q5(),
        D2Q9(),
        # D2Q17(),
    ]

    quadrature = last(quadratures)


    # NOTE we should base τ on ν and the speed of sound for the given qquadrature
    @show tgv
    N = tgv.NX
    τ = quadrature.speed_of_sound_squared * tgv.ν + 0.5
    # τ = 0.5;

    # initialize
    f_out = initialize!(quadrature, tgv, N)
    f_in = copy(f_out)

    stats = DataFrame(
        [Float64[], Float64[], Float64[], Float64[], Float64[]],
        [:density, :momentum, :total_energy, :kinetic_energy, :internal_energy]
    )
    process!(quadrature, f_in, 0, stats)

    # return f_in
    @inbounds for t = 0 : 2N
        if mod(t, 10) == 0
            @show t, t / 2N
        end
        process!(quadrature, f_in, t, stats)

        f_in = stream(quadrature, f_out)

        f_out = collide(SRT(τ), quadrature, f_in)
    end

    @show stats


    plot(
        plot(stats.density),
        plot(stats.momentum),
        plot(stats.total_energy),
        plot(stats.kinetic_energy),
        plot(stats.internal_energy)
    )

    gui()

    f_in, stats
end

# N = 2^5
example = TaylorGreenVortexExample(
    1.0 / 6.0,
    1
)

@time result = siumlate(example)

using Test

quadratures = [
    D2Q4(),
    D2Q5(),
    D2Q9()
]
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
