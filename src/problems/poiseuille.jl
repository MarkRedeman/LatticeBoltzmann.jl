struct PoiseuilleFlow{T <: Real, Int <: Integer} <: FluidFlowProblem
    rho_0::T
    u_max::T
    ν::T
    NX::Int
    NY::Int
    k::T
    domain_size::Tuple{T,T}
    G::T
end

function PoiseuilleFlow(
    ν = 1.0 / 6.0,
    scale = 2,
    NX = 5 * scale + 0,
    NY = NX,
    domain_size = (1.0, 1.0);
    static = true,
)
    # Since our velocity profile is constant over x it is sufficient to only take
    # 3 nodes in the x direction
    NX = 3
    u_max = 0.1 / scale

    return PoiseuilleFlow(1.0, u_max, ν, NX, NY, 1.0, domain_size, 1.0)
end

function density(
    q::Quadrature,
    problem::PoiseuilleFlow,
    x::T,
    y::T,
    timestep::Real = 0.0,
) where { T <: Real }
    return 1.0
end

function pressure(
    q::Quadrature,
    problem::PoiseuilleFlow,
    x::T,
    y::T,
    timestep::Real = 0.0,
) where { T <: Real }
    return 1.0
end

function velocity(problem::PoiseuilleFlow, x::T, y::T, timestep::Real = 0.0) where { T <: Real}
    return [
        y * (problem.domain_size[2] - y) * (problem.G / 2)
        0.0
    ]
end
function force(problem::PoiseuilleFlow, x::Int, y::Int, time::Real = 0.0) where { Int <: Integer }
    return [
        viscosity(problem) * problem.G
        0.0
    ]
end

function delta_x(problem::PoiseuilleFlow)
    return problem.domain_size[2] / problem.NY
end

has_external_force(problem::PoiseuilleFlow) = true

"""
Apply a bottom and top wall bounce back boundary condition
"""
boundary_conditions(problem::PoiseuilleFlow) = [
    BounceBack(North(), 1:problem.NX, 1:problem.NY),
    BounceBack(South(), 1:problem.NX, 1:problem.NY),
]
