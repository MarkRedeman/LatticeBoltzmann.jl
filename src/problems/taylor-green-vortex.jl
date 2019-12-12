struct TaylorGreenVortex <: lbm.FluidFlowProblem
    rho_0::Float64
    u_max::Float64
    ν::Float64
    NX::Int64
    NY::Int64
    domain_size::Tuple{Float64,Float64}
    static::Bool
    A::Float64
    B::Float64
    a::Float64
    b::Float64
end
function TaylorGreenVortex(
    ν = 1.0 / 6.0,
    scale = 2,
    NX = 16 * scale,
    NY = NX,
    domain_size = (2pi, 2pi);
    static = true,
    A = 1, B = -1, a = 1, b = 1
)
    u_max = sqrt(0.001) / scale
    u_max = 0.01 / scale
    Re = NX * u_max / ν
    @show Re

    p = TaylorGreenVortex(1.0, u_max, ν, NX, NY, domain_size, static, A, B, a, b)
    @show p

    if (p.u_max > sqrt(2 / 3) * delta_x(p) / delta_t(p))
        @warn p.u_max, sqrt(2 / 3) * delta_x(p) / delta_t(p)
    end
    p
end

function density(
    q::Quadrature,
    problem::TaylorGreenVortex,
    x::Float64,
    y::Float64,
    timestep::Float64 = 0.0,
)
    return pressure(q, problem, x, y, timestep)

    # If not athermal
    return 1.0
end

function pressure(
    q::Quadrature,
    problem::TaylorGreenVortex,
    x::Float64,
    y::Float64,
    timestep::Float64 = 0.0,
)
    # return 1.0
    a = problem.a
    A = problem.A
    b = problem.b
    B = problem.B

    P = -(1 / 4) *
        problem.rho_0 *
        decay(problem, x, y, timestep)^2 *
        (A^2 * cos(2 * a * x) + B^2 * cos(2 * b * y))

    return 1.0 + q.speed_of_sound_squared * problem.u_max^2 * P
end

function velocity(
    problem::TaylorGreenVortex,
    x::Float64,
    y::Float64,
    timestep::Float64 = 0.0,
)
    a = problem.a
    A = problem.A
    b = problem.b
    B = problem.B
    return decay(problem, x, y, timestep) * [
        A * cos(a * x) * sin(b * y)
        B * sin(a * x) * cos(b * y)
    ]
end

function velocity_gradient(
    problem::TaylorGreenVortex,
    x::Float64,
    y::Float64,
    timestep::Float64 = 0.0,
)
    a = problem.a
    A = problem.A
    b = problem.b
    B = problem.B

    u_x = -a * A * sin(x) * sin(y)
    v_y = -b * B * sin(x) * sin(y)
    u_y = b * A * cos(x) * cos(y)
    v_x = a * B * cos(x) * cos(y)

    return decay(problem, x, y, timestep) * [u_x v_x; u_y v_y]
end

function decay(problem::TaylorGreenVortex, x::Float64, y::Float64, timestep::Float64)
    a = problem.a
    b = problem.b

    return problem.static ?
        1.0 :
        exp(-(a^2 + b^2) * viscosity(problem) * timestep)
end

function force(problem::TaylorGreenVortex, x::Float64, y::Float64, time::Float64 = 0.0)
    return problem.static ?
        2 * viscosity(problem) * velocity(problem, x, y, 0.0) :
        [0.0 0.0]
end

has_external_force(problem::TaylorGreenVortex) = problem.static

# t_c(problem) = 1
# t_c(problem::TaylorGreenVortex) = ln(2) / (problem.ν * (problem.k_x^2 + problem.k_y^2))
