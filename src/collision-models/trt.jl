struct TRT <: CollisionModel
    τ_symmetric
    τ_asymmetric
end
# Magic parameter
function TRT(τ_symmetric, Δt, ν)
    # \Nabla = (τ_symmetric / Δt - .5) * (τ_asymmetric / Δt - .5)
    # should equal 1 / 4
end
