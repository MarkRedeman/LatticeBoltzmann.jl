"""
3. Veloctiy and stress only (ρ_0 = 1.0)
f₀(x) = f_eq(ρ₀(x), u₀(x), p₀(x))
ρ₀ = 1.0
u₀ = u₀(x)
p₀ = unkown
"""
struct ConstantDensity <: InitializationStrategy end

function initial_condition(::ConstantDensity, q::Quadrature, problem::FluidFlowProblem, x::T, y::T) where { T <: Real }
    u = lattice_velocity(q, problem, x, y)

    return equilibrium(q, one(T), u, one(T))
end
