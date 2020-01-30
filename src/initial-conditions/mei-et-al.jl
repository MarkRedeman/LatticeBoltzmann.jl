"""
4. Initialization scheme from Mei et al
"""
struct IterativeInitializationMeiEtAl{T <: Real} <: InitializationStrategy
    τ::T
    ϵ::T
end
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

struct DensityConvergence{FT <: Real, MT <: AbstractMatrix{FT}} <: StopCriteria
    ϵ::FT
    ρ_old::MT
    ρ::MT
end
function should_stop!(stop_criteria::DensityConvergence, q, f_in)
    nx, ny = size(stop_criteria.ρ_old)
    for x_idx = nx, y_idx = ny
        stop_criteria.ρ_old[x_idx, y_idx] = stop_criteria.ρ[x_idx, y_idx]
        stop_criteria.ρ[x_idx, y_idx] = density(q, f_in[x_idx, y_idx, :])
    end

    δρ = norm(stop_criteria.ρ - stop_criteria.ρ_old)

    # Stop procesing if the density converged or the method seems to be diverging
    return δρ < stop_criteria.ϵ || δρ > 100.0
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
        @info "Stopping due to convergence"
        return true
    end

    return false
end


struct ItirativeInitializationCollisionModel{ T <: Real, UT <: AbstractArray{T, 3}} <: CollisionModel
    τ::T

    # The incompressible equilibrium funciton can be split in a linear and a nonlinear component
    # the linear component depends on the local density while the nonlinear component depends
    # on the local velocity.
    # Since this collision operator keeps a constant velocity for each lattice node we can
    # precompute the nonlinear term for some performance improvements
    nonlinear_term::UT
end
function ItirativeInitializationCollisionModel(q::Quadrature, τ::T,  problem) where { T <: Real }
    nx, ny, nf = problem.NX, problem.NY, length(q.weights)
    nonlinear_term = Array{T}(undef, nx, ny, nf)
    ρ_0 = one(T)

    # Initialize f with the given velocity field
    x_range, y_range = range(problem)
    for x_idx = 1:nx, y_idx = 1:ny
        u = lattice_velocity(q, problem, x_range[x_idx], y_range[y_idx])
        u_squared = u[1]^2 + u[2]^2

        for f_idx = 1:nf
            u_dot_xi = dot(q.abscissae[:, f_idx], u)

            cs = q.speed_of_sound_squared
            a_H_1 = cs * u_dot_xi
            a_H_2 = cs^2 * (u_dot_xi * u_dot_xi)  - cs * u_squared

            nonlinear_term[x_idx, y_idx, f_idx] = q.weights[f_idx] * ρ_0 * (a_H_1 + a_H_2 / 2)
        end
    end

    return ItirativeInitializationCollisionModel(τ, nonlinear_term)
end
function collide!(
    collision_model::ItirativeInitializationCollisionModel,
    q::Quadrature,
    f_in,
    f_out;
    time = 0.0,
)
    τ = collision_model.τ
    nx, ny, nf = size(f_out)
    for x_idx = 1:nx, y_idx = 1:ny
        ρ = density(q, f_in[x_idx, y_idx, :])
        @inbounds for f_idx = 1:nf
            feq = q.weights[f_idx] * ρ + collision_model.nonlinear_term[x_idx, y_idx, f_idx]

            f_out[x_idx, y_idx, f_idx] = (1 - 1 / τ) * f_in[x_idx, y_idx, f_idx] + (1 / τ) * feq
        end
    end
    return
end
