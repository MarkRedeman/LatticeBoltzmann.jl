struct TGV <: lbm.FluidFlowProblem
    q::Quadrature
    ρ_0::Float64
    u_max::Float64
    u_0::Float64
    τ::Float64
    ν::Float64
    NX::Int64
    NY::Int64
    static::Bool
    domain_size::Tuple{Float64,Float64}
end
function TGV(
    q::Quadrature,
    τ,
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
    x::Float64,
    y::Float64,
    t::Float64 = 0.0,
)
    k_x = 2pi / problem.NX
    k_y = 2pi / problem.NY

    x *= problem.NX
    y *= problem.NY

    ν = (problem.τ - 0.5) / (problem.q.speed_of_sound_squared)
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
    x::Float64,
    y::Float64,
    t::Float64 = 0.0,
)
    k_x = 2pi / problem.NX
    k_y = 2pi / problem.NY

    x *= problem.NX
    y *= problem.NY

    ν = (problem.τ - 0.5) / (problem.q.speed_of_sound_squared)
    t_d = 1 / (ν * (k_x^2 + k_y^2))

    p_0 = 1.0
    ρ = 1.0


    return p_0 - ρ * (problem.u_0^2 / 4) * (
        (k_y / k_x) * cos(2 * k_x * x) +
        (k_x / k_y) * cos(2 * k_y * y)
    ) * exp(-2 * t / t_d)
end

function velocity(
    problem::TGV,
    x::Float64,
    y::Float64,
    t::Float64 = 0.0,
)
    k_x = 2pi / problem.NX
    k_y = 2pi / problem.NY

    x *= problem.NX
    y *= problem.NY

    ν = (problem.τ - 0.5) / (problem.q.speed_of_sound_squared)
    t_d = 1 / (ν * (k_x^2 + k_y^2))

    return problem.u_0 * exp(-t / t_d) * [
        - sqrt(k_y / k_x) * cos(k_x * x) * sin(k_y * y)
        sqrt(k_x / k_y) * sin(k_x * x) * cos(k_y * y)
    ]
end

function velocity_gradient(
    problem::TGV,
    x::Float64,
    y::Float64,
    t::Float64 = 0.0,
)
    k_x = 2pi / problem.NX
    k_y = 2pi / problem.NY

    x *= problem.NX
    y *= problem.NY

    ν = (problem.τ - 0.5) / (problem.q.speed_of_sound_squared)
    t_d = 1 / (ν * (k_x^2 + k_y^2))

    u_x = sqrt(k_y * k_x) * sin(k_x * x) * sin(k_y * y)
    v_y = - sqrt(k_y * k_x) * sin(k_x * x) * sin(k_y * y)

    u_y = - sqrt(k_y^3 / k_x) * cos(k_x * x) * cos(k_y * y)
    v_x = sqrt(k_x^3 / k_y) * cos(k_x * x) * cos(k_y * y)

    return exp(-t / t_d) * problem.u_0 * [u_x v_x; u_y v_y]
end

function decay(problem::TGV, x::Float64, y::Float64, t::Float64)
    ν = (problem.τ - 0.5) / (problem.q.speed_of_sound_squared)
    t_d = 1 / (ν * (k_x^2 + k_y^2))

    return problem.static ? 1.0 : exp(-t / t_d)
end

function force(problem::TGV, x::Float64, y::Float64, t::Float64 = 0.0)
    ν = (problem.τ - 0.5) / (problem.q.speed_of_sound_squared)

    return problem.static ? 2 * ν * velocity(problem, x, y, 0.0) : [0.0 0.0]
end

has_external_force(problem::TGV) = problem.static

viscosity(problem::TGV) = (problem.τ - 0.5) / (problem.q.speed_of_sound_squared)
delta_x(problem::TGV) = 1.0
delta_t(problem::TGV) = 1.0
