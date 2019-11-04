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
    p = DecayingShearFlow(1.0, u_max, ν_lb, NX, NY, 1.0, domain_size, static)
    p2 = DecayingShearFlow(1.0, u_max * scale, ν_lb * scale, NX, NY, 1.0, domain_size, static)

    # @show delta_t(p), delta_x(p), viscosity(p)
    # @show delta_t(p2), delta_x(p2), viscosity(p2)
    # return p2
    return DecayingShearFlow(1.0, u_max, ν_lb, NX, 3, 1.0, domain_size, static)
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

    return exp(-1.0 * problem.k^2 * viscosity(problem) * time )
end

function force(problem::DecayingShearFlow, x::Float64, y::Float64, time::Float64 = 0.0)
    if ! problem.static
        return [0.0 0.0]
    end

    A = 1.0
    B = 1.0

    return [
        0.0
        viscosity(problem) * problem.k^2 * B * cos(problem.k * (x - A * time))
    ]
end

has_external_force(problem::DecayingShearFlow) = problem.static
