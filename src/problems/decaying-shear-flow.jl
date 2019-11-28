struct DecayingShearFlow <: lbm.FluidFlowProblem
    rho_0::Float64
    u_max::Float64
    ν::Float64
    NX::Int64
    NY::Int64
    k::Float64
    domain_size::Tuple{Float64, Float64}
    static::Bool
    A::Float64
    B::Float64
end
function DecayingShearFlow(
    ν_lb = 1.0 / 6.0 , scale = 2, NX = 16 * scale, NY = NX, domain_size = (2pi, 2pi);
    static = true,
    A = 0.0,
    B = 1.0,
)
    u_max = 0.0025 / scale
    u_max = 0.02 / scale
    # ν_lb = ν_lb * scale
    Re = NX * u_max / ν_lb
    @show Re

   
    p = DecayingShearFlow(1.0, u_max, ν_lb, NX, 3, 1.0, domain_size, static, A, B)
    @show (delta_t(p), 1/delta_t(p))
    @show (delta_x(p), 1/delta_x(p))
    @show (p.u_max, 1/p.u_max)
    @show (p.domain_size[1] / delta_x(p)) / 16
    return p
end

function density(q::Quadrature, problem::DecayingShearFlow, x::Float64, y::Float64, time::Float64 = 0.0)
    return 1.0
    return pressure(q, problem, x, y, time)
end

function pressure(q::Quadrature, problem::DecayingShearFlow, x::Float64, y::Float64, time::Float64 = 0.0)
    return 1.0
    A = problem.A
    B = problem.B

    p = 1 + problem.B * (0.025) * q.speed_of_sound_squared * problem.u_max^2 * B * sin(problem.k * (x - A * time))^2 * decay(problem, x, y, time)^2
    return p
end
function pressure_tensor(q::Quadrature, problem::DecayingShearFlow, x::Float64, y::Float64, time::Float64 = 0.0)
    A = problem.A
    B = problem.B

    p = pressure(q, problem, x, y, time) * [1.0 0.0; 0.0 1.0]
    return p + deviatoric_tensor(q, problem, x, y, time)
    return p - q.speed_of_sound_squared * problem.k * B * sin(problem.k * x - problem.k * A * time) *
        decay(problem, x, y, time) * [
            0.0 1.0
            1.0 0.0
        ]
end
function deviatoric_tensor(q::Quadrature, problem::DecayingShearFlow, x::Float64, y::Float64, time::Float64 = 0.0)
    a = acceleration(problem, x, y, time)

    ν = problem.ν
    # ν = viscosity(problem)

    return ν * [
        2 * a[1, 1] a[1, 2] + a[2, 1]
        a[1, 2] + a[2, 1] 2 * a[2, 2]
    ]
end

function velocity(problem::DecayingShearFlow, x::Float64, y::Float64, time::Float64 = 0.0)
    A = problem.A
    B = problem.B

    k_y = 0.0
    k_x = problem.k

    return [
        A * cos(k_y * y - k_y * B * time) * decay(problem, x, y, time)
        B * cos(k_x * x - k_x * A * time) * decay(problem, x, y, time)
    ]
end
function acceleration(problem::DecayingShearFlow, x::Float64, y::Float64, time::Float64 = 0.0)
    A = problem.A
    B = problem.B

    k_y = 0.0
    k_x = problem.k

    u_x = 0.0
    u_y = A * k_y * sin(k_y * (y - B * time))
    v_x = B * k_x * sin(k_x * (x - A * time))
    v_y = 0.0

    return decay(problem, x, y, time) * [u_x v_x; u_y v_y]
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

    A = problem.A
    B = problem.B

    k_y = 0.0
    k_x = problem.k

    return [
        0.0
        viscosity(problem) * k_x^2 * B * cos(k_x * x - k_x * A * time)
    ]
end

has_external_force(problem::DecayingShearFlow) = problem.static
