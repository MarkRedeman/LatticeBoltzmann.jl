module Example
module TaylorGreenVortex

using Plots

struct TaylorGreenVortexExample end

include("../quadratures/quadrature.jl")
include("../quadratures/D2Q5.jl")
include("../quadratures/D2Q9.jl")
include("../stream.jl")
include("../collision.jl")

density(x, y) = 1.0
velocity(x, y) = [
    cos(x) * sin(y),
    sin(x) * cos(y)
]
function process!(quadrature, f_in, t, stats)
    if (mod(t, 1) == 0)
        # Density
        f_ρ = density(quadrature, f_in)

        # Momentum
        jx, jy = momentum(quadrature, f_in)


        # Velocity componetns
        u = jx ./ f_ρ
        v = jy ./ f_ρ

        s = (1200, 1200)
        velocity_field = contour(u[:, :, 1].^2 .+ v[:, :, 1].^2, fill=true, cbar=true, size=s, title="Momentum")
        N = size(u, 1)
        X = [i for i in range(1, size(u, 1), step = 1), j in range(1, size(u, 2), step = 1)]
        Y = [j for i in range(1, size(u, 1), step = 1), j in range(1, size(u, 2), step = 1)]
        # @show "process: ", u, v,
        @show "Conserved? ", sum(f_ρ), sum(jx + jy), sum(jx.^2 .+ jy.^2)

        quiver!(
            velocity_field,
            X, Y,
            # quiver=(u[:, :, 1], v[:, :, 1])
            quiver=(x, y) -> (u[x] / sqrt(u[x]^2 + v[x]^2), v[x] / sqrt(u[x]^2 + v[x]^2)),
            color="white",
            arrow=arrow(0.1, 0.1)
        )

        # @show u.^2 .+ v.^2
        # @show u
        plot(
            contour(f_ρ[:, :, 1], fill=true, clims=(0, 1.05), cbar=true, size=s),
            velocity_field,
            contour(u[:, :, 1].^2, fill=true, cbar=true, size=s, title="u"),
            contour(v[:, :, 1].^2, fill=true, cbar=true, size=s, title="v"),
            # size=(2 * 900, 600)
        )
        gui()
    end
end

function initialize!(quadrature, f_out)
    N = size(f_out, 1)
    # initial fields
    density_field = [
        density(x, y) for x in range(0, 2pi, length = N), y in range(0, 2pi, length = N)
    ]

    velocity_field = [
        0.2 * velocity(x, y) for x in range(0, 2pi, length = N), y in range(0, 2pi, length = N)
        # [0.2 0.2] for x in range(0, 2pi, length = N), y in range(0, 2pi, length = N)
        # 0.05 * rand(2) for x in range(0, 2pi, length = N), y in range(0, 2pi, length = N)
    ]

    # set initial conditions
    for x_idx = 1 : size(f_out, 1), y_idx = 1 : size(f_out, 2)
        # f_in[x_idx, y_idx, :] = equilibrium(D2Q9(), ρ, u, T)
        f_out[x_idx, y_idx, :] = equilibrium(
            quadrature,
            density_field[x_idx, y_idx],
            velocity_field[x_idx, y_idx],
            1.0 # Fow now we use a constant temperature
        )
    end

    jx, jy = momentum(quadrature, f_out)
    @show jx
    @show jy
end

function siumlate(::TaylorGreenVortexExample;)
    τ = 1.0;

    N = 2^4 + 1
    grid_size = (N, N)
    q = 9

    # initialize
    f_out = fill(0.0, grid_size..., 9)
    f_in = copy(f_out)

    # Idea: introduce an Initial Value Problem
    # problem = InitialValueProblem
    # solution = solve(problem, LBM(Lattice, CollisionModel))
    # LBM(Lattice, CollisionModel) can be a solution method


    stats = []

    quadrature = D2Q9()
    initialize!(quadrature, f_out)
    f_in = f_out
    process!(quadrature, f_in, 0, stats)

    # return f_in
    @inbounds for t = 1 : 10N
        process!(quadrature, f_in, t, stats)

        f_in = stream(quadrature, f_out)

        # collide
        # collide(lattice, collision_model)
        f_out = collide(SRT(τ), quadrature, f_in)
    end

    f_in
end

example = TaylorGreenVortexExample()

@time result = siumlate(example)

end
end


# Some thoughs: hide the storage of the distributions f inside of an interface
# so that we can do:
# for x in xs, y in ys, z in zs
#  f = fs[x, y ,z]::Vector
#
# the interface could hide it by returning
# return @view f_internal[x, y, z, :]
