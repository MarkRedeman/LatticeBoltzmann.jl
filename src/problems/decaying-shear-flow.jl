struct DecayingShearFlow <: lbm.FluidFlowProblem
    rho_0::Float64
    u_max::Float64
    ν::Float64
    NX::Int64
    NY::Int64
    domain_size::Tuple{Float64, Float64}
    static::Bool
    A::Float64
    B::Float64
    k_x::Float64
    k_y::Float64
end
function DecayingShearFlow(
    ν_lb = 1.0 / 6.0 , scale = 2, NX = 8 * scale, NY = NX, domain_size = (2pi, 2pi);
    static = true,
    A = 1.0,
    B = 1.0,
    k_x = 1.0,
    k_y = 0.0
)
    u_max = 0.02 / scale
    Re = NX * u_max / ν_lb
    @show Re

    if (k_y == 0.0)
        NY = 3
    end
    if (k_x == 0.0)
        NX = 3
    end

    return DecayingShearFlow(1.0, u_max, ν_lb, NX, NY, domain_size, static, A, B, k_x, k_y)
end

function density(q::Quadrature, problem::DecayingShearFlow, x::Float64, y::Float64, time::Float64 = 0.0)
    return 1.0
    return pressure(q, problem, x, y, time)
end

function pressure(q::Quadrature, problem::DecayingShearFlow, x::Float64, y::Float64, time::Float64 = 0.0)
    return 1.0
    A = problem.A
    B = problem.B

    k_y = problem.k_y
    k_x = problem.k_x

    p = 1 + problem.B * (0.025) * q.speed_of_sound_squared * problem.u_max^2 * B * sin(k_x * (x - A * time))^2 * decay(problem, x, y, time)^2
    return p
end

function velocity(problem::DecayingShearFlow, x::Float64, y::Float64, time::Float64 = 0.0)
    A = problem.A
    B = problem.B

    k_y = problem.k_y
    k_x = problem.k_x
    if problem.static
        return [
            A * cos(k_y * y - k_y * B * time)
            B * cos(k_x * x - k_x * A * time)
        ]
    end

    return [
        A * cos(k_y * y - k_y * B * time) * exp(-1.0 * k_y^2 * viscosity(problem) * time )
        B * cos(k_x * x - k_x * A * time) * exp(-1.0 * k_x^2 * viscosity(problem) * time )
    ]
end
function velocity_gradient(problem::DecayingShearFlow, x::Float64, y::Float64, time::Float64 = 0.0)
    A = problem.A
    B = problem.B

    k_y = problem.k_y
    k_x = problem.k_x

    u_x = 0.0
    u_y = A * k_y * sin(k_y * (y - B * time)) * exp(-1.0 * k_y^2 * viscosity(problem) * time )
    v_x = B * k_x * sin(k_x * (x - A * time)) * exp(-1.0 * k_x^2 * viscosity(problem) * time )
    v_y = 0.0

    if problem.static
        u_y = A * k_y * sin(k_y * (y - B * time))
        v_x = B * k_x * sin(k_x * (x - A * time))
    end

    return [u_x v_x; u_y v_y]
end

function decay(problem::DecayingShearFlow, x::Float64, y::Float64, time::Float64)
    if (problem.static)
        return 1.0
    end

    A = problem.A
    B = problem.B

    k_y = problem.k_y
    k_x = problem.k_x

    return exp(-1.0 * k_x^2 * viscosity(problem) * time )
end

function force(problem::DecayingShearFlow, x::Float64, y::Float64, time::Float64 = 0.0)
    if ! problem.static
        return [0.0 0.0]
    end

    A = problem.A
    B = problem.B

    k_y = problem.k_y
    k_x = problem.k_x

    return [
        viscosity(problem) * k_y^2 * A * cos(k_y * y - k_y * B * time),
        viscosity(problem) * k_x^2 * B * cos(k_x * x - k_x * A * time)
    ]
end

has_external_force(problem::DecayingShearFlow) = problem.static
