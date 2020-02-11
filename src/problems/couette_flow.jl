struct CouetteFlow{T <: Real, Int <: Integer} <: FluidFlowProblem
    rho_0::T
    u_max::T
    ν::T
    NX::Int
    NY::Int
    domain_size::Tuple{T, T}
end

function CouetteFlow(
    ν = 1.0 / 6.0,
    scale = 2,
    NX = 5 * scale + 0,
    NY = NX,
    domain_size = (1.0, 1.0),
)
    # # Since our velocity profile is constant over x it is sufficient to only take
    # # 3 nodes in the x direction
    NX = 1
    u_max = 0.015 / scale
    u_max = 0.01 / scale
    # u_max = sqrt(0.001) / scale

    return CouetteFlow(1.0, u_max, ν, NX, NY, domain_size)
end

function density(
    q::Quadrature,
    problem::CouetteFlow,
    x::T,
    y::T,
    timestep::Real = 0.0,
) where {T <: Real}
    return 1.0
end

function pressure(
    q::Quadrature,
    problem::CouetteFlow,
    x::T,
    y::T,
    timestep::Real = 0.0,
) where {T <: Real}
    return 1.0
end

function velocity(problem::CouetteFlow, x::T, y::T, timestep::Real = 0.0) where {T <: Real}
    return [
        y
        0.0
    ]
end
function velocity_gradient(
    problem::CouetteFlow,
    x::T,
    y::T,
    timestep::Real = 0.0,
) where {T <: Real}
    u_x = 0.0
    v_y = 0.0
    u_y = 1.0
    v_x = 0.0

    return [u_x v_x; u_y v_y]
end
function force(
    problem::CouetteFlow,
    x::Int,
    y::Int,
    time::Real = 0.0,
) where {Int <: Integer}
    return [
        0.0
        0.0
    ]
end

"""
Apply a bounce back to the bottom of the domain and a moving wall to the
top of the domain
"""
boundary_conditions(problem::CouetteFlow) = [
    BounceBack(South(), 1:(problem.NX), 1:(problem.NY)),
    MovingWall(North(), 1:(problem.NX), 1:(problem.NY), [problem.u_max, 0]),
]
