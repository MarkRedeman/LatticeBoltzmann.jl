export TaylorGreenVortexDecay,
    density,
    velocity,
    temperature,
    decay,
    force

struct TaylorGreenVortexDecay <: lbm.InitialValueProblem
    scale
    rho_0
    u_max
    ν
    NX
    NY
    k_x
    k_y

    function TaylorGreenVortexDecay(ν = 1.0 / 6.0 , scale = 2, NX = 16 * scale, NY = NX)
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
function density(tgv::TaylorGreenVortexDecay, x::Int, y::Int, timestep::Int = 0)
    X = x + 0.5
    Y = y + 0.5
    u_max = tgv.u_max
    kx = tgv.k_x
    ky = tgv.k_y

    P = -0.25 * tgv.rho_0 * u_max * u_max * (
        (ky / kx) * cos(2.0 * kx * X) + (kx / ky)*cos(2.0 * ky * Y)
    ) * decay(tgv, x, y, timestep)^2;

    return 1.0 + 3.0 * P
end
function velocity(tgv::TaylorGreenVortexDecay, x::Int, y::Int, timestep::Int = 0)
    X = x + 0.5
    Y = y + 0.5
    u_max = tgv.u_max
    kx = tgv.k_x
    ky = tgv.k_y

    return decay(tgv, x, y, timestep) .* [
        -u_max * sqrt(ky / kx) * cos(kx * X) * sin(ky * Y),
        u_max * sqrt(kx / ky) * sin(kx * X) * cos(ky * Y)
    ]
end
function decay(tgv::TaylorGreenVortexDecay, x, y, timestep::Int)
    td = 1.0 / (tgv.ν * (tgv.k_x^2 + tgv.k_y^2));

    exp(-1.0 * timestep / td)
end

function temperature(tgv::TaylorGreenVortexDecay, x::Int, y::Int, timestep::Int = 0)
    1.0
end

force(t::TaylorGreenVortexDecay, x, y) = 0
force(t, x, y) = (t.a^2 + t.b^2) * (t.length^2 / t.Re) * velocity(t, x, y)
