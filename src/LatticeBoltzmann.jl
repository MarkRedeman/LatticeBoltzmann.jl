module LatticeBoltzmann

import Base: range

using LinearAlgebra
using Plots

include("quadratures.jl")
include("boundary_conditions.jl")
include("problems/problems.jl")
include("initial_conditions.jl")
include("stream.jl")
include("collision_models.jl")
include("analysis/stopping_criteria.jl")
include("analysis/convergence.jl")
include("process.jl")
include("initial_conditions/mei_et_al.jl")

include("interop/plots.jl")

export CollisionModel, SRT, TRT, MRT
export Quadrature, D2Q4, D2Q5, D2Q9, D2Q13, D2Q17, D2Q21, D2Q37, opposite

export analayze_convergence
# Problems
export DecayingShearFlow
export LidDrivenCavityFlow
export TaylorGreenVortex
export CouetteFlow
export GenericFluidFlowProblem
export LinearizedThermalDiffusion, LinearizedTransverseShearWave
export PoiseuilleFlow

export process!,
    apply_boundary_conditions!,
    density,
    velocity,
    pressure,
    temperature,
    decay,
    force,
    FluidFlowProblem,
    viscosity,
    delta_t
export Lattice,
    initialize,
    AnalyticalEquilibriumAndOffEquilibrium,
    AnalyticalEquilibrium,
    AnalyticalVelocity,
    IterativeInitialization,
    # Thermodynamics
    density,
    momentum,
    pressure,
    total_energy,
    kinetic_energy,
    internal_energy,
    temperature,
    dimension,
    # Equilibria
    equilibrium,
    equilibrium!,
    hermite_equilibrium,
    hermite_first_nonequilibrium,
    hermite,
    # Siumlation
    stream!,
    collide!,
    simulate

struct LatticeBoltzmannMethod{Q<:Quadrature,T,CM<:CollisionModel,PM<:ProcessingMethod,BCs <: AbstractVector{<:BoundaryCondition}}
    f_stream::T
    f_collision::T
    quadrature::Q
    collision_model::CM
    boundary_conditions::BCs
    processing_method::PM
end
function LatticeBoltzmannMethod(
    problem,
    quadrature;
    collision_model = SRT,
    initialization_strategy = InitializationStrategy(problem),
    process_method
)
    f_stream = initialize(initialization_strategy, quadrature, problem, collision_model)
    f_collision = copy(f_stream)

    LatticeBoltzmannMethod(
        f_stream,
        f_collision,
        quadrature,
        CollisionModel(collision_model, quadrature, problem),
        boundary_conditions(problem),
        process_method,
    )
end
function simulate(
    problem::FluidFlowProblem,
    q::Quadrature;
    process_method = nothing,
    should_process = true,
    initialization_strategy = InitializationStrategy(problem),
    t_end = 1.0,
    collision_model = SRT,
)
    Δt = delta_t(problem)
    n_steps = round(Int, t_end / Δt)

    if isnothing(process_method)
        process_method = ProcessingMethod(problem, should_process, n_steps)
    end

    model = LatticeBoltzmannMethod(
        problem,
        q,
        collision_model = collision_model,
        initialization_strategy = initialization_strategy,
        process_method = process_method,
    )

    simulate(model, 0:n_steps)
end
function simulate(model::LatticeBoltzmannMethod, time)
    Δt = isdefined(model.processing_method, :problem) ? delta_t(model.processing_method.problem) : 0.0

    @inbounds for t in time
        collide!(model, time = t * Δt)
        stream!(model)
        apply_boundary_conditions!(model, time = t * Δt)

        if process_step!(model, t + 1)
            return model
        end
    end

    process_step!(model, last(time) + 1)

    model
end

# The following are helper functions which use the old interface where we
# explicitely pass the quadrature, collision model etc.
# This will likely be refactored so that we can write specialized function
# for each specific LatticeBoltzmannModel (once we also introduce ddf models)

function collide!(model::LatticeBoltzmannMethod; time)
    collide!(
        model.collision_model,
        model.quadrature,
        f_new = model.f_collision,
        f_old = model.f_stream,
        time = time,
    )
end

function stream!(model::LatticeBoltzmannMethod)
    stream!(model.quadrature, f_new = model.f_stream, f_old = model.f_collision)
end

function apply_boundary_conditions!(model::LatticeBoltzmannMethod; time = 0.0)
    apply!(
        model.boundary_conditions,
        model.quadrature,
        model.f_stream,
        model.f_collision,
        time = time,
    )
end

function process_step!(model::LatticeBoltzmannMethod, t::Int64)
    next!(model.processing_method, model.quadrature, model.f_stream, t)
end

end
