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

function force(tgv::TaylorGreenVortexDecay, x::Int, y::Int, time::Int = 0)
    @show "HOIHOI"
    return (1.0 / (tgv.ν * (tgv.k_x^2 + tgv.k_y^2))) * [0.01, 0.4] # velocity(tgv, x, y, 0)
    # return (t.a^2 + t.b^2) * (t.length^2 / t.Re) * velocity(tgv, x, y, 0)
# (1.0 / (tgv.ν * (tgv.k_x^2 + tgv.k_y^2))) * t.length / (t.Re / t.speed)

#     Decay: - (t.a^2 + t.b^2) * (t.speed * t.length) / (t.Re)
#     - 1.0 / (tgv.ν * (tgv.k_x^2 + tgv.k_y^2));

    # Re = NX * u_max / ν
    # Re = t.speed * t.length / ν
    # td = - (t.a^2 + t.b^2) * ν
    #     -u_max * sqrt(ky / kx) * cos(kx * X) * sin(ky * Y),
    #     u_max * sqrt(kx / ky) * sin(kx * X) * cos(ky * Y)

    # td = - 1.0 / (tgv.ν * (tgv.k_x^2 + tgv.k_y^2));
    # td = - (t.a^2 + t.b^2) * (t.speed * t.length) / (t.Re)
    # * timestep

    # t.A = u_max * sqrt(ky / kx)
    # t.B = -u_max * sqrt(ky / kx)
    # t.a * t.length * x = kx * X
    # t.b * t.length * y = ky * Y
        # -u_max * sqrt(ky / kx) * cos(kx * X) * sin(ky * Y),
        # u_max * sqrt(kx / ky) * sin(kx * X) * cos(ky * Y)
#         t.A * cos(t.a * t.length * x)sin(t.b * t.length * y),
#         t.B * sin(t.a * t.length * x)cos(t.b * t.length * y)
    return (t.a^2 + t.b^2) * (t.length^2 / t.Re) * velocity(tgv, x, y, 0)
end
# force(t, x, y) = (t.a^2 + t.b^2) * (t.length^2 / t.Re) * velocity(t, x, y)

# function velocity(t::TaylorGreenVortex, x, y, time)
#     return decay(t, time * t.length / t.speed) * (1 / t.speed) * [
#         t.A * cos(t.a * t.length * x)sin(t.b * t.length * y),
#         t.B * sin(t.a * t.length * x)cos(t.b * t.length * y)
#     ]
# end
# decay(t::TaylorGreenVortex, time) = exp(- (t.a^2 + t.b^2) * (t.speed * t.length) / (t.Re) * time)
# decay(t::StaticVortex, time) = 1
# decay(t::DecayingVortex, time) = exp(- (t.a^2 + t.b^2) * (t.speed * t.length) / (t.Re) * time)
