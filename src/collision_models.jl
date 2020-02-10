abstract type CollisionModel end

include("collision_models/srt.jl")
include("collision_models/trt.jl")
include("collision_models/mrt.jl")

"""
Our default collision model uses the Single Relaxation Time method
"""
function CollisionModel(
    cm::Type{<:CollisionModel},
    q::Quadrature,
    problem::FluidFlowProblem,
)
    return CollisionModel(SRT, q, problem)
end
CollisionModel(cm::CollisionModel, ::Quadrature, ::FluidFlowProblem) = cm
collide!(cm, q::Quadrature; time, f_new, f_old) = collide!(cm, q, f_old, f_new, time = time)
