"""
4. Initialization scheme from Mei et al
"""
struct IterativeInitializationMeiEtAl{T <: Real} <: InitializationStrategy
    τ::T
    ϵ::T
end
IterativeInitializationMeiEtAl() = IterativeInitialization()
IterativeInitialization() = IterativeInitializationMeiEtAl(1.0, 1E-7)

function initialize(
    strategy::IterativeInitializationMeiEtAl,
    q::Quadrature,
    problem::FluidFlowProblem,
    cm = SRT;
    process_method = ProcessingMethod(problem, true, 1000)
)
    # Initialize f such that it has a constant density for all nodes, We don't
    # care about the velocity here since this will be set in the collision method
    f = [
        q.weights[f_idx]
        for x_idx = 1:problem.NX, y_idx = 1:problem.NY, f_idx = 1:length(q.weights)
    ]

    # Simultae the lattice botlzmann method with a special collision operator that uses
    # u_0(x) as the velocity of each lattice node
    # The method stops if the density has converged
    model = LatticeBoltzmannMethod(
        copy(f),
        copy(f),
        q,
        ItirativeInitializationCollisionModel(q, strategy.τ, problem),
        boundary_conditions(problem),
        ProcessIterativeInitialization(strategy.ϵ, problem, process_method)
    )

    simulate(model, 1:10000)

    return model.f_stream
end

struct ProcessIterativeInitialization3{T <: ProcessingMethod} <: ProcessingMethod
    stop_criteria::StopCriteria
    internal_process_method::T
    n_steps::Int
end
function ProcessIterativeInitialization(ϵ, problem, process_method)
    nx, ny = problem.NX, problem.NY

    ρ = zeros(nx, ny)
    ρ_old = zeros(nx, ny)

    ProcessIterativeInitialization3(
        DensityConvergence(ϵ, ρ, ρ_old),
        process_method,
        100
    )
end
function next!(process_method::ProcessIterativeInitialization3{T}, q, f_in, t) where {T}
    if mod(t, max(10, round(Int, process_method.n_steps / 25))) == 0
        # next!(process_method.internal_process_method, q, f_in, 0)
    end

    if (should_stop!(process_method.stop_criteria, q, f_in))
        return true
    end

    return false
end
