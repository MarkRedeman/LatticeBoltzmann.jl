
abstract type StopCriteria end
struct NoStoppingCriteria <: StopCriteria end
mutable struct MeanVelocityStoppingCriteria3 <: StopCriteria
    old_mean_velocity::Float64
    tolerance::Float64
    problem::InitialValueProblem
end
StopCriteria(problem::InitialValueProblem) = NoStoppingCriteria()
StopCriteria(problem::PoiseuilleFlow) = MeanVelocityStoppingCriteria3(0.0, 1e-12, problem)

should_stop!(::StopCriteria, q, f_in) = false
function should_stop!(stop_criteria::MeanVelocityStoppingCriteria3, q, f_in)
    f = Array{Float64}(undef, size(f_in, 3))
    u = zeros(dimension(q))
    Nx = size(f_in, 1)
    Ny = size(f_in, 2)

    divide_by = 0
    u_mean = 0.0

    @inbounds for x_idx in 1:Nx, y_idx in 2:Ny-1
        if ! is_fluid(stop_criteria.problem, x_idx, y_idx)
            continue
        end
        divide_by += 1

        # Calculated
        @inbounds for f_idx = 1 : size(f_in, 3)
            f[f_idx] = f_in[x_idx, y_idx, f_idx]
        end
        ρ = density(q, f)
        velocity!(q, f, ρ, u)

        u_mean += u[1]
    end

    u_mean /= divide_by

    converged = abs(u_mean / stop_criteria.old_mean_velocity - 1)

    if (converged < stop_criteria.tolerance)
        return true
    end

    stop_criteria.old_mean_velocity = u_mean

    return false
end
