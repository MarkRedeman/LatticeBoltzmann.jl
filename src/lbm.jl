module lbm

import Base: range

using LinearAlgebra
using Plots

include("quadratures.jl")
include("boundary-conditions.jl")
include("problems/problems.jl")
include("initial-conditions.jl")
include("stream.jl")
include("collision.jl")
include("analysis.jl")
include("process.jl")
include("initial-conditions/mei-et-al.jl")

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

    lbm = LatticeBoltzmannMethod(
        problem,
        q,
        collision_model = collision_model,
        initialization_strategy = initialization_strategy,
        process_method = process_method,
    )

    simulate(lbm, 0:n_steps)
end
function simulate(lbm::LatticeBoltzmannMethod, time)
    Δt = isdefined(lbm.processing_method, :problem) ? delta_t(lbm.processing_method.problem) : 0.0

    @inbounds for t in time
        collide!(lbm, time = t * Δt)
        stream!(lbm)
        apply_boundary_conditions!(lbm, time = t * Δt)

        if process_step!(lbm, t + 1)
            break
        end
    end

    process_step!(lbm, last(time) + 1)

    lbm
end

# The following are helper functions which use the old interface where we
# explicitely pass the quadrature, collision model etc.
# This will likely be refactored so that we can write specialized function
# for each specific LatticeBoltzmannModel (once we also introduce ddf models)

function collide!(lbm::LatticeBoltzmannMethod; time)
    collide!(
        lbm.collision_model,
        lbm.quadrature,
        f_new = lbm.f_collision,
        f_old = lbm.f_stream,
        time = time,
    )
end

function stream!(lbm::LatticeBoltzmannMethod)
    stream!(lbm.quadrature, f_new = lbm.f_stream, f_old = lbm.f_collision)
end

function apply_boundary_conditions!(lbm::LatticeBoltzmannMethod; time = 0.0)
    apply!(
        lbm.boundary_conditions,
        lbm.quadrature,
        lbm.f_stream,
        lbm.f_collision,
        time = time,
    )
end

function process_step!(lbm::LatticeBoltzmannMethod, t::Int64)
    next!(lbm.processing_method, lbm.quadrature, lbm.f_stream, t)
end

end
