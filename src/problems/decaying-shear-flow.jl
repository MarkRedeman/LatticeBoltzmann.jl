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
    ν = 1.0 / 6.0 , scale = 2, NX = 16 * scale, NY = NX, domain_size = (2pi, 2pi); static = true
)
    u_max = 0.02 / scale
    @show u_max
    Re = NX * u_max / ν
    @show Re
    return DecayingShearFlow(1.0, u_max, ν, NX, NY, 1.0, domain_size, static)
end

function viscosity(problem::DecayingShearFlow)
    return problem.ν * (2pi)^2 * (1 / problem.NX^2 + 1 / problem.NY^2) * 0.5
    return problem.ν
end

function density(q::Quadrature, problem::DecayingShearFlow, x::Float64, y::Float64, timestep::Float64 = 0.0)
    return 1.0
end

function pressure(q::Quadrature, problem::DecayingShearFlow, x::Float64, y::Float64, timestep::Float64 = 0.0)
    return 1.0
end

function velocity(problem::DecayingShearFlow, x::Float64, y::Float64, timestep::Float64 = 0.0)
    u_max = problem.u_max
    v_max = problem.u_max

    return [
       v_max
       u_max * cos(problem.k * (x - v_max * timestep)) * decay(problem, x, y, timestep)
    ]
end
function decay(problem::DecayingShearFlow, x::Float64, y::Float64, timestep::Float64)
    if (problem.static)
        return 1.0
    end
    ν = problem.ν * (1 / problem.NX^2 + 0 / problem.NY^2)
    ν = problem.ν * (2pi)^2 * (1 / problem.NX^2 + 1 / problem.NY^2) * 0.5
    ν = viscosity(problem)
    Δt = delta_t(problem)

    return exp(-1.0 * problem.k^2 * ν * timestep / Δt)
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
    u_max = problem.u_max
    v_max = problem.u_max

    ν = problem.ν * (1 / problem.NX^2 + 0 / problem.NY^2)
    ν = problem.ν * (2pi)^2 * (1 / problem.NX^2 + 1 / problem.NY^2) * 0.5
    ν = viscosity(problem)
    Δt = delta_t(problem)

    decay = exp(-1.0 * problem.k^2 * ν * time / Δt)
    decay = 1.0

    # δt = delta_t(problem)
    δt = 0.5 * (2pi)^2 * (1 / problem.NX^2 + 1 / problem.NY^2)

    return δt * problem.ν  * [
        0.0
        u_max  * problem.k^2 * cos(problem.k * x - problem.k * v_max * time)
    ] * decay
    return 0.0 * velocity(problem, x, y, time)
end

function delta_t(problem::DecayingShearFlow)
    # Note: is actually δx * u_max
    return Δt = (2pi)^1 * (1 / problem.NX + 1 / problem.NY) * 0.5
    return 1.0
    ν = viscosity(problem)
    return 1.0
    return ν
    Δt = (2pi)^2 * (1 / problem.NX^2 + 1 / problem.NY^2)# * (problem.k_x^2 + problem.k_y^2)

    return Δt
end

function delta_x(problem::DecayingShearFlow)
    return (2pi) * (1 / problem.NX + 1 / problem.NY)
    return (2pi)^2 * (1 / problem.NX^2 + 0 / problem.NY^2)
end

function lattice_viscosity(problem::InitialValueProblem)
    Re = 1.0 / problem.ν
    δ_t = 1.0
    δ_x = problem.domain_size[1] / problem.NX
    δ_t / (δ_x^2 * Re)
end

has_external_force(problem::DecayingShearFlow) = problem.static
