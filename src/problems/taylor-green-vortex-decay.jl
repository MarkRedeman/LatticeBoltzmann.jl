export TaylorGreenVortexDecay,
    TaylorGreenVortexExample,
    density,
    velocity,
    pressure,
    temperature,
    decay,
    force,
    initialize

struct TaylorGreenVortexExample <: lbm.InitialValueProblem
    scale
    rho_0
    u_max
    ν
    NX
    NY
    k_x
    k_y

    function TaylorGreenVortexExample(ν = 1.0 / 6.0 , scale = 2, NX = 16 * scale, NY = NX)
        u_max = 0.02 / scale
        Re = NX * u_max / ν
        @show Re
        return new(
            scale,
            1.0,
            u_max,
            ν,
            NX,
            NY,
            2pi / NX,
            2pi / NY,
        )
    end
end

function initialize!(quadrature::Quadrature, tgv::TaylorGreenVortexExample)
    force_field = Array{Float64}(undef, tgv.NX, tgv.NY, dimension(quadrature))
    f = Array{Float64}(undef, tgv.NX, tgv.NY, length(quadrature.weights))

    # NOTE: we have periodic boundaries
    x_range = range(0, 2pi, length=tgv.NX + 1)
    y_range = range(0, 2pi, length=tgv.NY + 1)
    for x_idx in 1:tgv.NX, y_idx in 1:tgv.NY
        x = x_range[x_idx]
        y = y_range[y_idx]

        f[x_idx, y_idx, :] = equilibrium(
            quadrature,
            density(quadrature, tgv, x, y),
            velocity(tgv, x, y),
            pressure(quadrature, tgv, x, y) / density(quadrature, tgv, x, y)
        )
        force_field[x_idx, y_idx, :] = force(tgv, x, y)
    end

    τ = quadrature.speed_of_sound_squared * tgv.ν + 0.5
    collision_operator = SRT_Force(τ, force_field)
    collision_operator = SRT(τ)

    return f, collision_operator
end

function density(q::Quadrature, tgv::TaylorGreenVortexExample, x::Float64, y::Float64, timestep::Float64 = 0.0)
    return pressure(q, tgv, x, y, timestep)
    
    # If not athermal
    return 1.0
end

function pressure(q::Quadrature, tgv::TaylorGreenVortexExample, x::Float64, y::Float64, timestep::Float64 = 0.0)
    P = -0.25 * tgv.rho_0 * tgv.u_max^2 * (
        (tgv.k_y / tgv.k_x) * cos(2.0 * x) +
        (tgv.k_x / tgv.k_y) * cos(2.0 * y)
    ) * decay(tgv, x, y, timestep)^2;

    return 1.0 + q.speed_of_sound_squared * P
end

function velocity(tgv::TaylorGreenVortexExample, x::Float64, y::Float64, timestep::Float64 = 0.0)
    u_max = tgv.u_max

    return decay(tgv, x, y, timestep) .* [
      -u_max * sqrt(tgv.k_y / tgv.k_x) * cos(x) * sin(y),
       u_max * sqrt(tgv.k_x / tgv.k_y) * sin(x) * cos(y)
    ]
end
function decay(tgv::TaylorGreenVortexExample, x::Float64, y::Float64, timestep::Float64)
    return exp(-1.0 * timestep)
end

function force(tgv::TaylorGreenVortexExample, x::Float64, y::Float64, time::Float64 = 0.0)
    return (1 / tgv.ν) * (tgv.k_x^2 + tgv.k_y^2) * velocity(tgv, x, y, time)
end
