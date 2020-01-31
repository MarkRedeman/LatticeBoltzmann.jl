struct TGV{T <: Real, Int <: Integer} <: FluidFlowProblem
    q::Quadrature
    ρ_0::T
    u_max::T
    u_0::T
    τ::T
    ν::T
    NX::Int
    NY::Int
    static::Bool
    domain_size::Tuple{T,T}
end
function TGV(
    q::Quadrature,
    τ, # Effective relaxation time
    scale = 2,
    NX = 16 * scale,
    NY = NX,
    # u_max = sqrt(0.001) / scale
    u_max = 0.02 / scale
)
    ν = (τ - 0.5) / (q.speed_of_sound_squared)
    Re = NX * u_max / ν
    @show Re

    TGV(q, 1.0, 1.0, u_max, τ, ν, NX, NY, false, (1.0, 1.0))
end

function density(
    q::Quadrature,
    problem::TGV,
    x::T,
    y::T,
    t::Real = 0.0,
) where { T <: Real }
    k_x = 2pi / problem.NX
    k_y = 2pi / problem.NY

    x *= problem.NX
    y *= problem.NY

    ν = problem.ν
    t_d = 1 / (ν * (k_x^2 + k_y^2))

    return problem.ρ_0  * (
        1.0 - q.speed_of_sound_squared * (problem.u_0^2 / 4) * (
            (k_y / k_x) * cos(2 * k_x * x) +
            (k_x / k_y) * cos(2 * k_y * y)
        ) * exp(-2 * t / t_d)
    )
end

function pressure(
    q::Quadrature,
    problem::TGV,
    x::T,
    y::T,
    t::Real = 0.0,
) where { T <: Real }
    k_x = 2pi / problem.NX
    k_y = 2pi / problem.NY

    x *= problem.NX
    y *= problem.NY

    ν = problem.ν
    t_d = 1 / (ν * (k_x^2 + k_y^2))

    p_0 = problem.ρ_0
    ρ = 1.0

    return p_0 - q.speed_of_sound_squared * ρ * (problem.u_0^2 / 4) * (
        (k_y / k_x) * cos(2 * k_x * x) +
        (k_x / k_y) * cos(2 * k_y * y)
    ) * exp(-2 * t / t_d)
end

function velocity(
    problem::TGV,
    x::T,
    y::T,
    t::Real = 0.0,
) where { T <: Real }
    k_x = 2pi / problem.NX
    k_y = 2pi / problem.NY

    x *= problem.NX
    y *= problem.NY

    ν = problem.ν
    t_d = 1 / (ν * (k_x^2 + k_y^2))

    return problem.u_0 * exp(-t / t_d) * [
        - sqrt(k_y / k_x) * cos(k_x * x) * sin(k_y * y)
        sqrt(k_x / k_y) * sin(k_x * x) * cos(k_y * y)
    ]
end

function velocity_gradient(
    problem::TGV,
    x::T,
    y::T,
    t::Real = 0.0,
) where { T <: Real }
    k_x = 2pi / problem.NX
    k_y = 2pi / problem.NY

    x *= problem.NX
    y *= problem.NY

    ν = problem.ν
    t_d = 1 / (ν * (k_x^2 + k_y^2))

    u_x = sqrt(k_y * k_x) * sin(k_x * x) * sin(k_y * y)
    v_y = - sqrt(k_y * k_x) * sin(k_x * x) * sin(k_y * y)

    u_y = - sqrt(k_y^3 / k_x) * cos(k_x * x) * cos(k_y * y)
    v_x = sqrt(k_x^3 / k_y) * cos(k_x * x) * cos(k_y * y)

    return exp(-t / t_d) * problem.u_0 * [u_x v_x; u_y v_y]
end

function decay(problem::TGV, x::T, y::T, t::Real) where { T <: Real }
    ν = problem.ν
    k_x = 2pi / problem.NX
    k_y = 2pi / problem.NY
    t_d = 1 / (ν * (k_x^2 + k_y^2))

    return problem.static ? 1.0 : exp(-t / t_d)
end

function force(problem::TGV, x::T, y::T, t::Real = 0.0) where { T <: Real }
    ν = problem.ν

    return problem.static ? 2 * ν * velocity(problem, x, y, 0.0) : [0.0 0.0]
end

has_external_force(problem::TGV) = problem.static

viscosity(problem::TGV) = problem.ν # (problem.τ - 0.5) / (problem.q.speed_of_sound_squared)
delta_x(problem::TGV) = 1.0
delta_t(problem::TGV) = 1.0

function decay_time(problem)
    ν = viscosity(problem)
    k_x = 2pi / problem.NX
    k_y = 2pi / problem.NY
    t_d = 1 / (ν * (k_x^2 + k_y^2))
end
