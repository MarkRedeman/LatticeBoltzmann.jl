struct ProcessIterativeInitialization{T <: ProcessingMethod} <: ProcessingMethod
    stop_criteria::StopCriteria
    internal_process_method::T
    n_steps::Int
end
function ProcessIterativeInitialization(ϵ, problem, process_method)
    nx, ny = problem.NX, problem.NY

    ρ = zeros(nx, ny)
    ρ_old = zeros(nx, ny)

    ProcessIterativeInitialization(
        DensityConvergence(ϵ, ρ, ρ_old),
        process_method,
        100
    )
end
function next!(process_method::ProcessIterativeInitialization{T}, q, f_in, t) where {T}
    if mod(t, max(10, round(Int, process_method.n_steps / 25))) == 0
        # next!(process_method.internal_process_method, q, f_in, 0)
    end

    if (should_stop!(process_method.stop_criteria, q, f_in))
        return true
    end

    return false
end
