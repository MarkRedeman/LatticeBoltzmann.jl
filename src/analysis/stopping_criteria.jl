abstract type StopCriteria end
struct NoStoppingCriteria <: StopCriteria end
mutable struct MeanVelocityStoppingCriteria{ T <: Real } <: StopCriteria
    old_mean_velocity::T
    tolerance::T
    problem::FluidFlowProblem
end
StopCriteria(problem::FluidFlowProblem) = NoStoppingCriteria()
StopCriteria(problem::PoiseuilleFlow) = MeanVelocityStoppingCriteria(0.0, 1e-12, problem)
StopCriteria(problem::CouetteFlow) = MeanVelocityStoppingCriteria(0.0, 1e-7, problem)
StopCriteria(problem::LidDrivenCavityFlow) =
    MeanVelocityStoppingCriteria(0.0, 1e-5, problem)
StopCriteria(problem::DecayingShearFlow) =
    problem.static ? MeanVelocityStoppingCriteria(0.0, 1e-8, problem) : NoStoppingCriteria()

should_stop!(::StopCriteria, q, f_in) = false
function should_stop!(stop_criteria::MeanVelocityStoppingCriteria{T}, q, f_in) where { T }
    f = Array{T}(undef, size(f_in, 3))
    u = zeros(dimension(q))
    Nx = size(f_in, 1)
    Ny = size(f_in, 2)

    divide_by = 0
    u_mean = 0.0

    @inbounds for x_idx = 1:Nx, y_idx = 1:Ny
        divide_by += 1

        # Calculated
        @inbounds for f_idx = 1:size(f_in, 3)
            f[f_idx] = f_in[x_idx, y_idx, f_idx]
        end
        ρ = density(q, f)
        velocity!(q, f, ρ, u)

        u_mean += u[1]
    end

    u_mean /= divide_by

    converged = abs(u_mean / stop_criteria.old_mean_velocity - 1)

    if (converged < stop_criteria.tolerance)
        @info "Stopping due to convergence"
        return true
    end

    if (isnan(u_mean))
        @warn "nan in velocity profile"
        return true
    end

    stop_criteria.old_mean_velocity = u_mean

    return false
end

struct VelocityConvergenceStoppingCriteria{T <: Real, VT <: AbstractVector{T}} <: StopCriteria
    old_velocity::Array{VT, 2}
    tolerance::T
    problem::FluidFlowProblem
end
function VelocityConvergenceStoppingCriteria(tolerance::T, problem) where { T <: Real }
    return VelocityConvergenceStoppingCriteria(
        [zeros(T, 2) for x_idx = 1:problem.NX, y_idx = 1:problem.NY],
        tolerance,
        problem
    )
end

function should_stop!(stop_criteria::VelocityConvergenceStoppingCriteria{T}, q, f_in) where { T }
    f = Array{T}(undef, size(f_in, 3))
    u = zeros(dimension(q))
    Nx = size(f_in, 1)
    Ny = size(f_in, 2)

    old_velocity_norm = 0.0 # norm(stop_criteria.old_velocity)
    error = 0.0

    @inbounds for x_idx = 1:Nx, y_idx = 1:Ny
        # Calculated
        @inbounds for f_idx = 1:size(f_in, 3)
            f[f_idx] = f_in[x_idx, y_idx, f_idx]
        end

        ρ = density(q, f)
        velocity!(q, f, ρ, u)

        u_old = stop_criteria.old_velocity[x_idx, y_idx]

        error += (
            ((u[1] - u_old[1])^2 + (u[2] - u_old[2])^2)
        )

        old_velocity_norm += u_old[1]^2 + u_old[2]^2
        stop_criteria.old_velocity[x_idx, y_idx] = copy(u)
    end

    converged = sqrt(error) / old_velocity_norm

    if (converged < stop_criteria.tolerance)
        @info "Stopping due to convergence"
        return true
    end

    if (isnan(converged))
        @warn "nan in velocity profile"
        return true
    end

    # stop_criteria.old_mean_velocity = u_mean

    return false
end
