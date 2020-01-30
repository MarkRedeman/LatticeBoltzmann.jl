"""
3. Veloctiy and stress only (ρ_0 = 1.0)
f₀(x) = f_eq(ρ₀(x), u₀(x), p₀(x))
ρ₀ = 1.0
u₀ = u₀(x)
p₀ = p₀(x)
"""
struct AnalyticalVelocity <: InitializationStrategy end

# TODO compute the density
average_density(problem::FluidFlowProblem) = 1.0
average_pressure(problem::FluidFlowProblem) = 1.0
function density(
    q::Quadrature,
    problem::FluidFlowProblem,
    x,
    y,
    p̄ = average_pressure(problem),
    ρ̄ = average_density(problem)
)
    δρ = (pressure(problem, x, y) - p̄) / q.speed_of_sound_squared

    return ρ̄ + δρ
end
# TODO Move computing the averages to the initialize function so that these can be computed globaly
function initial_condition(::AnalyticalVelocity, q::Quadrature, problem::FluidFlowProblem, x::T, y::T) where { T <: Real }
    p̄ = average_pressure(problem)
    ρ̄ = average_density(problem)

    # Compute the density from the pressure, average pressure and average density
    ρ = density(q, problem, x, y, p̄, ρ̄)

    u = lattice_velocity(q, problem, x, y),
    temperature = pressure(problem, x, y) / ρ

    return equilibrium(q, ρ, u, temperature)
end
