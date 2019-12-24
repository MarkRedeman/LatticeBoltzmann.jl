"""
2. Velocity and pressure only (initialize with equilibrium)
f₀(x) = f_eq(ρ₀(x), u₀(x), p₀(x))
ρ₀ = ρ₀(x)
u₀ = u₀(x)
p₀ = p₀(x)
"""
struct AnalyticalEquilibrium <: InitializationStrategy end
function initial_condition(::AnalyticalEquilibrium, q::Quadrature, problem::FluidFlowProblem, x::Float64, y::Float64)
    ρ = lattice_density(q, problem, x, y)
    u = lattice_velocity(q, problem, x, y)
    T = lattice_temperature(q, problem, x, y)

    return equilibrium(q, ρ, u, T)
end
