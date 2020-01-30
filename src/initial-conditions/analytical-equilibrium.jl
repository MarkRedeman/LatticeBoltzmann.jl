"""
2. Velocity and pressure only (initialize with equilibrium)
f₀(x) = f_eq(ρ₀(x), u₀(x), p₀(x))
ρ₀ = ρ₀(x)
u₀ = u₀(x)
p₀ = p₀(x)
"""
struct AnalyticalEquilibrium <: InitializationStrategy end
function initial_condition(::AnalyticalEquilibrium, q::Quadrature, problem::FluidFlowProblem, x::T, y::T) where { T <: Real }
    return equilibrium(q, problem, x, y)
end
