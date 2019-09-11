module TaylorGreen

"""
Taylor Green Vortex in 2 dimensions
"""
abstract type TaylorGreenVortex{T}

struct DecayingVortex{T<:Real} <: TaylorGreenVortex{T}
    a::T
    b::T
    A::T
    B::T
    Re::T # Reynolds number
    length::T
    speed::T

    function DecayingVortex(a, b, A, B, ν, length, speed)
        a * A + b * B != 0 ? throw(ArgumentError("The given vortex is not incompressible")) : new(a, b, A, B, ν, length, speed)
    end
end

struct StaticVortex{T<:Real} <: TaylorGreenVortex{T}
    a::T
    b::T
    A::T
    B::T
    Re::T # Reynolds number
    length::T
    speed::T

    function StaticVortex(a, b, A, B, ν, length, speed)
        a * A + b * B != 0 ? throw(ArgumentError("The given vortex is not incompressible")) : new(a, b, A, B, ν, length, speed)
    end
end

# Constructors
TaylorGreenVortex{T<:Real}(a::T, b::T, A::T, B::T, ν::T, length::T, speed::T) = DecayingVortex{T}(a, b, A, B, ν, length, speed)
TaylorGreenVortex{T<:Real}(a::T, b::T, A::T, B::T, ν::T) = DecayingVortex{T}(a, b, A, B, ν, 1., 1.)
TaylorGreenVortex() = TaylorGreenVortex{Float64}(1., 1., 1., -1., 1.)

velocity(t::TaylorGreenVortex, x, y) = velocity(t, x, y, 0.0)
function velocity(t::TaylorGreenVortex, x, y, time)
    return decay(t, time * t.length / t.speed) * (1 / t.speed) * [
        t.A * cos(t.a * t.length * x)sin(t.b * t.length * y),
        t.B * sin(t.a * t.length * x)cos(t.b * t.length * y)
    ]
end

decay(t::TaylorGreenVortex, time) = exp(- (t.a^2 + t.b^2) * (t.speed * t.length) / (t.Re) * time)
decay(t::StaticVortex, time) = 1
decay(t::DecayingVortex, time) = exp(- (t.a^2 + t.b^2) * (t.speed * t.length) / (t.Re) * time)

pressure(t::TaylorGreenVortex, x, y) = pressure(t, x, y, 0.0)
function pressure(t::TaylorGreenVortex, x, y, time)
    return (t.A^2 / (4 * t.speed^2)) * (cos(2 * t.a * t.length * x) + (t.a^2 / t.b^2) * cos(2 * t.b * t.length * y)) * decay(t, time * t.length / t.speed)
end

force(t::DecayingVortex, x, y) = 0
force(t::StaticVortex, x, y) = (t.a^2 + t.b^2) * (t.length^2 / t.Re) * velocity(t, x, y)


end
