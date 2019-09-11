# http://www.timm-krueger.de/downloads/Krueger_Edmonton_scaling.pdf
# http://s467657437.online.de/wp-content/uploads/2019/08/Krueger_Edmonton_scaling.pdf
# http://s467657437.online.de/wp-content/uploads/2019/08/Krueger_Edmonton_IBM.pdf

module TGV_Test_Case

using JLD

include("taylor-green.jl")
include("lbm.jl")

abstract Simulation

immutable SycnSimulation <: Simulation
    problem
    lattice
    boundary_conditions
end

immutable ParallelSimulation <: Simulation
    problem
    lattice
    boundary_conditions
    processors
    current_processor
    messages #? worked in c++ each message would contain info on where to send the data (data, node_idx, f_idx)
end

# typealias LBM_Problem{Model} LBM.Problem{Model<:TaylorGreen.TaylorGreenVortex}
immutable LBM_Problem{Model<:TaylorGreen.TaylorGreenVortex}
    model::Model
    # model::Model
    # quadrature::Quadrature
    # collision::Collision
    N::Int64
    N_iter::Int64
end

include("plot-lbm.jl")

function initial_condition(problem::LBM_Problem, x, y)::LBM.Distribution
    Δx = 1. / problem.N
    Δt = 1. / problem.N_iter

    return LBM.equilibrium(
        1.0,
        (Δt / Δx) *
        # (Δt / Δx) * # when adding force it goes to this speed
        TaylorGreen.velocity(problem.model, x, y)
    )
end

# Rename to initialize(problem)::Lattice ::Simulation?
# May also include logic for parallel programming?
function initial_condition(problem::LBM_Problem)::Array{LBM.Distribution, 2}
    return [
        initial_condition(problem, x, y) for x = linspace(0.0, 1.0, problem.N), y = linspace(0.0, 1.0, problem.N)
    ]
end

function relaxation_rate(problem::LBM_Problem)
    const Δx = 1.0 / problem.N
    const Δt = 1.0 / problem.N_iter

    const ν = (Δt / Δx^2) * (1 / problem.model.Re)
    const ω = 1. / (3ν + 0.5)

    return ω
end

# Specific collision model with forcing
function collision_model{T}(problem::LBM_Problem{TaylorGreen.StaticVortex{T}})
    const Δx = 1.0 / problem.N
    const Δt = 1.0 / problem.N_iter

    force = (Δt^2 / Δx) * (Δt / Δx) * [
        TaylorGreen.force(problem.model, x, y) for x = linspace(0.0, 1.0, problem.N), y = linspace(0.0, 1.0, problem.N)
    ]

    const ω = relaxation_rate(problem)
    return (fs) -> LBM.collide.(fs, ω, force)
end

# General collision model
function collision_model(problem::LBM_Problem)
    const ω = relaxation_rate(problem)
    return (fs) -> LBM.collide.(fs, ω)
end

function simulate(problem::LBM_Problem, after_update)

    # GENERAL FLOW CONSTANTS
    drho = 0.001;
    delta_rho = -drho * (1 - 2.0 * rand(lx, ly));

    # INITIAL CONDITION FOR BOTH DISTRIBUTION FUNCTIONS: (T=0) ==> TIn(i) = t(i)
    fIn = zeros(9, lx, ly)
    gIn = zeros(9, lx, ly)
    for i = 1:9, x = 1:lx, y = 1:ly
        fIn[i, x, y] = tNS[i] .* (1.0 + delta_rho[x, y]);
        gIn[i, x, y] = tNS[i] .* (1.0 - delta_rho[x, y]);
    end

    # @sync
    fs = initial_condition(problem)
    collide = collision_model(problem)

    after_update(fs, 0)
    @inbounds for timestep = 1:problem.N_iter
        fs = collide(fs)

        # Check if this is faster...
        # fs .= collide.(fs)

        fs = LBM.stream(fs)

        # Send values from ghost nodes (must be after streaming since now some nodes might miss values)
        # @sync communicate()
        # apply_boundary_conditions() # streaming, bounce back etc.

        after_update(fs, timestep)
    end

    return fs;
end

# Note that if I set a = 2. then we will have twice as many voritces
# if the solution is unstable it will tend towards 4 vortices (4 yellow, 4 purple)
a = 1.
A = 1.
length = 2.0π
speed = A

cells = [4, 8, 16, 32, 64, 128, 256, 512]
cells = [32, 64, 128, 256, 512]
cells = [8, 16, 32, 64]

# It seems that taking a timestep such that u_lb < 0.03 gives stable results
timesteps = [
    (cells) -> 500
    (cells) -> 50 * cells
    (cells) -> 50 * cells^2
]

# for Re = [1e-1, 1e0, 1e1, 1e2, 1e3]
    Re = 2e2
    model = TaylorGreen.DecayingVortex{Float64}(a, a, A, -A, Re, length, speed)

    compute_convergence_rate(
        model,
        cells,
        (cells) -> 500 * cells
    )
# end
end
