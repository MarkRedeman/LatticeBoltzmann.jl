# http://www.timm-krueger.de/downloads/Krueger_Edmonton_scaling.pdf

module TGV_Test_Case

# using JLD

include("taylor-green.jl")
include("lbm.jl")

abstract type Simulation end

struct SycnSimulation <: Simulation
    problem
    lattice
    boundary_conditions
end

struct ParallelSimulation <: Simulation
    problem
    lattice
    boundary_conditions
    processors
    current_processor
    messages #? worked in c++ each message would contain info on where to send the data (data, node_idx, f_idx)
end

# typealias LBM_Problem{Model} LBM.Problem{Model<:TaylorGreen.TaylorGreenVortex}
struct LBM_Problem{Model<:TaylorGreen.TaylorGreenVortex}
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

function error(problem::LBM_Problem, fs, x, y, t = 1.0)
    Δx = 1. / problem.N
    Δt = 1. / problem.N_iter

    const grid = linspace(0.0, 1.0, problem.N)
    exact = (x, y) -> TaylorGreen.velocity(problem.model, grid[x], grid[y], t)
    u = LBM.velocity(fs[x, y])

    return norm((Δx / Δt) * u - exact(x, y))
end


function analyze!(problem, fs, errors)
    lx = problem.N
    steps = problem.N_iter

    # Keep track of the L_2 error and the error on a local point
    total_error = sqrt(sum([error(problem, fs, x, y)^2 for x = 1:lx, y = 1:lx] / (lx^2)))
    local_error = error(problem, fs, Int(ceil(lx / 2)), Int(ceil(lx / 2)))

    # Keep track of the error and the parameters used to get this error
    push!(errors, Dict(
        :reynolds    => problem.model.Re,
        :N           => problem.N,
        :N_iter      => problem.N_iter,
        :total_error => total_error,
        :local_error => local_error
    ))

    # Plot the convergence rate using the currently known error rates
    l, = size(errors)
    plot(
        [errors[i][:N] for i = 1:l],
        [
            [errors[i][:local_error] for i = 1:l],
            [4.0^(-i) for i = 1:l]
        ],
        yscale = :log10, xscale = :log10,
        label = ["experiment", "expected"],
        xaxis = "N",
        yaxis = "error",
        marker = (:square)
    )

    gui()
end


function plot_taylor_green_after_update(problem)
    stats = initial(PlotStatistics, problem)
    const tPlot = Int(ceil(problem.N_iter / 10))

    return function(fs, timestep)
        update!(stats, fs)
        if timestep % tPlot == 0
            plot_taylor(fs, timestep, stats)
        end

        # Check how much we've decayed
        if timestep == 1
            # Store results?

            # Show decay rate of the supposedly static solution
            if (typeof(problem.model) <: TaylorGreen.StaticVortex)
                final = last(stats.decay_rate.decay_rate_of_experiment.decay)
                initial = first(stats.decay_rate.decay_rate_of_experiment.decay)
                # @show initial, final, initial / final
            end
        end
    end
end

function store_solutions_in_dir(problem)
    # We will first create a special folder within the solutions folder for
    # this particular problem
    # Next we will store the problem and its intermediate solutions in that folder

    path = "solutions/reynolds_$(problem.model.Re)/$(problem.N)_$(problem.N_iter)"

    mkpath(path)
    save(string(path, "/problem.jld"), "problem", problem)

    save_after = cld(problem.N_iter, 10)

    return function(fs, timestep)
        if timestep == 0 || timestep == problem.N_iter # || timestep % save_after == 0
            save(string(path, "/solution_", timestep, ".jld"), "fs", fs)
        end
    end
end



function compute_convergence_rate(
    model,
    amount_of_cells::Array{Int, 1} = [4, 8, 16, 32, 64],
    steps = (cells) -> 25cells
)
    errors = []

    for cells ∈ amount_of_cells
        problem = LBM_Problem(model, cells, steps(cells))

        @time fs = simulate(
            problem,
            plot_taylor_green_after_update(problem)
            # store_solutions_in_dir(problem)
        );

        # Update our errors based on new solution
        analyze!(problem, fs, errors)
    end

    @show errors

    save("solutions/errors_reynolds_$(model.Re).jld", "errors", errors)

    return errors
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
