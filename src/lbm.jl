module lbm

using DataFrames
using BenchmarkTools
using TimerOutputs

include("quadratures.jl")
include("boundary-conditions.jl")
include("problems/problems.jl")
include("stream.jl")
include("collision.jl")
include("analysis.jl")
include("process.jl")

export Quadrature, Lattice, D2Q4, D2Q5, D2Q9, D2Q17,
    opposite,
    initialize,
    density,
    momentum,
    pressure,
    total_energy,
    kinetic_energy,
    internal_energy,
    temperature,
    dimension,
    equilibrium,
    equilibrium!,
    hermite_equilibrium,
    hermite_first_nonequilibrium,
    hermite,
    stream,
    stream!,
    CollisionModel,
    SRT,
    TRT,
    collide!,
    simulate

struct LatticeBoltzmannMethod{Q<:Quadrature, T, CM <: CollisionModel, PM <: ProcessingMethod}
    f_stream::T
    f_collision::T
    quadrature::Q
    collision_model::CM
    boundary_conditions::Vector{<:BoundaryCondition}
    processing_method::PM
end
function LatticeBoltzmannMethod(
    problem,
    quadrature,
    collision_model;
    n_steps = 100,
    should_process = false
)
    f_stream = initialize(quadrature, problem, collision_model)
    f_collision = similar(f_stream)

    processing_method =
    return LatticeBoltzmannMethod(
        f_stream,
        f_collision,
        quadrature,
        CollisionModel(
            collision_model,
            quadrature,
            problem
        ),
        boundary_conditions(
            problem
        ),
        ProcessingMethod(
            problem,
            should_process,
            n_steps
        )
    )
end

function siumlate(
    problem::FluidFlowProblem,
    q::Quadrature;
    process_method = nothing,
    should_process = true,
    t_end = 1.0,
    collision_model = SRT
)
    Δt = delta_t(problem)
    n_steps = round(Int, t_end / Δt)

    lbm = LatticeBoltzmannMethod(
        problem,
        q,
        collision_model,
        n_steps = n_steps,
        should_process = should_process
    )

    @inbounds for t = 0:n_steps
        if process_step!(lbm, t)
            break
        end

        collide!(lbm, time = t * Δt, problem = problem)
        stream!(lbm)
        apply_boundary_conditions!(lbm, time = t * Δt)
    end

    process_step!(lbm, n_steps)

    lbm
end

# The following are helper functions which use the old interface where we
# explicitely pass the quadrature, collision model etc.
# This will likely be refactored so that we can write specialized function
# for each specific LatticeBoltzmannModel (once we also introduce ddf models)

function collide!(lbm; time, problem)
    collide!(
        lbm.collision_model,
        lbm.quadrature,
        f_new = lbm.f_collision,
        f_old = lbm.f_stream,
        time = time,
        problem = problem
    )
end

function stream!(lbm)
    stream!(
        lbm.quadrature,
        f_new = lbm.f_stream,
        f_old = lbm.f_collision
    )
end

function apply_boundary_conditions!(lbm; time = 0.0)
    apply!(
        lbm.boundary_conditions,
        lbm.quadrature,
        lbm.f_stream,
        lbm.f_collision,
        time = time
    )
end

function process_step!(lbm, t::Int64)
    next!(lbm.processing_method, lbm.quadrature, lbm.f_stream, t)
end

end
