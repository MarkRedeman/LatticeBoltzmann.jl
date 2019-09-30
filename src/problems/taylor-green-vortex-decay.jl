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
    rho_0::Float64
    u_max::Float64
    ν::Float64
    NX::Int64
    NY::Int64
    k_x::Float64
    k_y::Float64
    domain_size::Tuple{Float64, Float64}
    static::Bool
end
function TaylorGreenVortexExample(ν = 1.0 / 6.0 , scale = 2, NX = 16 * scale, NY = NX, domain_size = (2pi, 2pi); static = true)
    u_max = 0.02 / scale
    Re = NX * u_max / ν
    @show Re
    return TaylorGreenVortexExample(
        1.0,
        u_max,
        ν,
        NX,
        NY,
        domain_size[1] / NX,
        domain_size[2] / NY,
        domain_size,
        static
    )
end

function viscosity(problem::TaylorGreenVortexExample)
    return problem.ν * (problem.k_x^2 + problem.k_y^2)
end

function density(q::Quadrature, tgv::TaylorGreenVortexExample, x::Float64, y::Float64, timestep::Float64 = 0.0)
    # return pressure(q, tgv, x, y, timestep)
    
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

    return decay(tgv, x, y, timestep) * [
      -u_max * sqrt(tgv.k_y / tgv.k_x) * cos(x) * sin(y),
       u_max * sqrt(tgv.k_x / tgv.k_y) * sin(x) * cos(y)
    ]
end
function decay(tgv::TaylorGreenVortexExample, x::Float64, y::Float64, timestep::Float64)
    return tgv.static ? 1.0 : exp(-1.0 * timestep)
end

function force(problem::TaylorGreenVortexExample, x_idx::Int64, y_idx::Int64, time::Float64 = 0.0)
    x_range = range(0, problem.domain_size[1], length=problem.NX + 1)
    y_range = range(0, problem.domain_size[2], length=problem.NY + 1)

    x = x_range[x_idx]
    y = y_range[y_idx]

    return force(problem, x, y, time)
end
function force(tgv::TaylorGreenVortexExample, x::Float64, y::Float64, time::Float64 = 0.0)
    # return tgv.u_max * (1 / tgv.ν) * (tgv.k_x + tgv.k_y)^2 * velocity(tgv, x, y, 0.0)
    return tgv.static ? delta_t(tgv) * velocity(tgv, x, y, 0.0) : [0.0 0.0]

    return tgv.ν * (tgv.k_x^2 + tgv.k_y^2)
end

function delta_t(problem::TaylorGreenVortexExample)
    ν = viscosity(problem)
    Δt = ν #* (problem.k_x^2 + problem.k_y^2)

    return Δt
end

has_external_force(problem::TaylorGreenVortexExample) = problem.static
