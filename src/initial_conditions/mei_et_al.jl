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
    model = LatticeBoltzmannModel(
        copy(f),
        copy(f),
        q,
        IterativeInitializationCollisionModel(q, strategy.τ, problem),
        boundary_conditions(problem),
        ProcessIterativeInitialization(strategy.ϵ, problem, process_method)
    )

    simulate(model, 1:10000)

    return model.f_stream
end

