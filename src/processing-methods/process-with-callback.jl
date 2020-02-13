struct ProcessWithCallback{Callback} <: ProcessingMethod
    problem::FluidFlowProblem
    every_t::Union{Int, Vector{Int}}
    callback::Callback
end

function next!(pm::TakeSnapshots, q, f_in, t::Int64)
    if pm.every_t isa Int && mod(t, pm.every_t) != 0
        return false
    end
    if pm.every_t isa Vector{Int} && t ∉ pm.every_t
        return false
    end

    pm.callback(q, f_in, t)

    return false
end
