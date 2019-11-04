export TaylorGreenVortexDecay, TaylorGreenVortexExample

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
    p = TaylorGreenVortexExample(1.0, u_max, ν, NX, NY, domain_size[1] / NX, domain_size[2] / NY, domain_size, static)

    if (p.u_max > sqrt(2/3) * delta_x(p) / delta_t(p))
        @warn p.u_max, sqrt(2/3) * delta_x(p) / delta_t(p)
    end
    p
end

function density(q::Quadrature, tgv::TaylorGreenVortexExample, x::Float64, y::Float64, timestep::Float64 = 0.0)
    # return pressure(q, tgv, x, y, timestep)
    
    # If not athermal
    return 1.0
end

function pressure(q::Quadrature, tgv::TaylorGreenVortexExample, x::Float64, y::Float64, timestep::Float64 = 0.0)
    # return 1.0
    P = -0.25 * tgv.rho_0 * (
        (tgv.k_y / tgv.k_x) * cos(2.0 * x) +
        (tgv.k_x / tgv.k_y) * cos(2.0 * y)
    ) * decay(tgv, x, y, timestep)^2

    return 1.0 + q.speed_of_sound_squared * tgv.u_max^2 * P
end

function velocity(tgv::TaylorGreenVortexExample, x::Float64, y::Float64, timestep::Float64 = 0.0)
    u_max = tgv.u_max
    u_max = 1.0

    # if (timestep == 0.0 && tgv.static)
    #     return [
    #         0.0
    #         0.0
    #     ]
    # end

    return decay(tgv, x, y, timestep) * [
      -u_max * sqrt(tgv.k_y / tgv.k_x) * cos(x) * sin(y),
       u_max * sqrt(tgv.k_x / tgv.k_y) * sin(x) * cos(y)
    ]
end
function decay(tgv::TaylorGreenVortexExample, x::Float64, y::Float64, timestep::Float64)
    return tgv.static ? 1.0 : exp(-2.0 * timestep * viscosity(tgv))
end

function force(tgv::TaylorGreenVortexExample, x::Float64, y::Float64, time::Float64 = 0.0)
    return tgv.static ? 2viscosity(tgv) * velocity(tgv, x, y, 0.0) : [0.0 0.0]
end

has_external_force(problem::TaylorGreenVortexExample) = problem.static
