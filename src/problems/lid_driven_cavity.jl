struct LidDrivenCavityFlow{T <: Real, Int <: Integer} <: FluidFlowProblem
    rho_0::T
    u_max::T
    ν::T
    NX::Int
    NY::Int
    domain_size::Tuple{T, T}
end

function LidDrivenCavityFlow(
    ν = 1.0 / 6.0,
    scale = 2,
    NX = 16 * scale,
    NY = NX,
    domain_size = (1.0, 1.0),
)
    # # Since our velocity profile is constant over x it is sufficient to only take
    # # 3 nodes in the x direction
    # NX = 5
    u_max = 0.01 / scale

    return LidDrivenCavityFlow(1.0, u_max, ν, NX, NY, domain_size)
end

function density(
    q::Quadrature,
    problem::LidDrivenCavityFlow,
    x::T,
    y::T,
    timestep::Real = 0.0,
) where {T <: Real}
    return 1.0
end

function pressure(
    q::Quadrature,
    problem::LidDrivenCavityFlow,
    x::T,
    y::T,
    timestep::Real = 0.0,
) where {T <: Real}
    return 1.0
end

function velocity(
    problem::LidDrivenCavityFlow,
    x::T,
    y::T,
    timestep::Real = 0.0,
) where {T <: Real}
    return [
        0.0
        0.0
    ]
end

"""
Apply a moving wall to the top of the domain and bounce back to all other
sides of the domain
"""
boundary_conditions(problem::LidDrivenCavityFlow) = [
    BounceBack(East(), 1:(problem.NX), 1:(problem.NY)),
    BounceBack(South(), 1:(problem.NX), 1:(problem.NY)),
    BounceBack(West(), 1:(problem.NX), 1:(problem.NY)),
    MovingWall(North(), 1:(problem.NX), 1:(problem.NY), [problem.u_max, 0]),
]
