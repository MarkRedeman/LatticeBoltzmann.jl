export DecayingShearFlow

struct DecayingShearFlow <: lbm.InitialValueProblem
    rho_0::Float64
    u_max::Float64
    ν::Float64
    NX::Int64
    NY::Int64
    k::Float64
    domain_size::Tuple{Float64, Float64}
    static::Bool
end
function DecayingShearFlow(
    ν_lb = 1.0 / 6.0 , scale = 2, NX = 16 * scale, NY = NX, domain_size = (2pi, 2pi); static = true
)
    u_max = 0.02 / scale
    # ν_lb = ν_lb * scale
    @show u_max
    Re = NX * u_max / ν_lb
    @show Re
    return DecayingShearFlow(1.0, u_max, ν_lb, NX, NY, 1.0, domain_size, static)
end

# Dimensionless
function viscosity(problem::DecayingShearFlow)
    return problem.ν * delta_x(problem)^2 / delta_t(problem)
end

function density(q::Quadrature, problem::DecayingShearFlow, x::Float64, y::Float64, timestep::Float64 = 0.0)
    return 1.0
end

function pressure(q::Quadrature, problem::DecayingShearFlow, x::Float64, y::Float64, timestep::Float64 = 0.0)
    return 1.0
end

function velocity(problem::DecayingShearFlow, x::Float64, y::Float64, time::Float64 = 0.0)
    A = 1.0
    B = 1.0

    return [
        A
        B * cos(problem.k * (x - A * time)) * decay(problem, x, y, time)
    ]
end
function decay(problem::DecayingShearFlow, x::Float64, y::Float64, time::Float64)
    if (problem.static)
        return 1.0
    end

    ν = viscosity(problem)
    return exp(-1.0 * problem.k^2 * ν * time )
end

function force(problem::DecayingShearFlow, x_idx::Int64, y_idx::Int64, time::Float64 = 0.0)
    x_range = range(0, problem.domain_size[1], length=problem.NX + 1)
    y_range = range(0, problem.domain_size[2], length=problem.NY + 1)

    x = x_range[x_idx]
    y = y_range[y_idx]

    return force(problem, x, y, time)
end
function force(problem::DecayingShearFlow, x::Float64, y::Float64, time::Float64 = 0.0)
    if ! problem.static
        return [0.0 0.0]
    end
    A = 1.0
    B = 1.0

    ν = viscosity(problem)
    return [
        0.0
        ν * problem.k^2 * B * cos(problem.k * x - problem.k * A * time)
    ]
end

function delta_t(problem::DecayingShearFlow)
    return delta_x(problem) * problem.u_max
end

function delta_x(problem::DecayingShearFlow)
    return problem.domain_size[1] * (1 / problem.NX)
end

function lattice_viscosity(problem::InitialValueProblem)
    return problem.ν
    Re = 1.0 / problem.ν
    δ_t = 1.0
    δ_x = problem.domain_size[1] / problem.NX
    δ_t / (δ_x^2 * Re)
end

has_external_force(problem::DecayingShearFlow) = problem.static

lattice_density(q, problem::DecayingShearFlow, x, y, t = 0.0) = density(q, problem, x, y, t)
lattice_velocity(q, problem::DecayingShearFlow, x, y, t = 0.0) = problem.u_max * velocity(problem, x, y, t)
lattice_pressure(q, problem::DecayingShearFlow, x, y, t = 0.0) = pressure(q, problem, x, y, t)
lattice_force(problem::DecayingShearFlow, x, y, t = 0.0) = problem.u_max * delta_t(problem) * force(problem, x, y, t)
