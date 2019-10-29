export PoiseuilleFlow

import Base: range

struct PoiseuilleFlow <: lbm.InitialValueProblem
    rho_0::Float64
    u_max::Float64
    ν::Float64
    NX::Int64
    NY::Int64
    k::Float64
    domain_size::Tuple{Float64, Float64}
end

function PoiseuilleFlow(
    ν = 1.0 / 6.0, scale = 2, NX = 5 * scale + 0, NY = NX, domain_size = (1.0, 1.0) ; static = true
)
    # Since our velocity profile is constant over x it is sufficient to only take
    # 3 nodes in the x direction
    NX = 3
    u_max = 0.1 / scale
   
    return PoiseuilleFlow(
        1.0,
        u_max,
        ν,
        NX,
        NY,
        1.0,
        domain_size
    )
end

function density(q::Quadrature, problem::PoiseuilleFlow, x::Float64, y::Float64, timestep::Float64 = 0.0)
    return 1.0
end

function pressure(q::Quadrature, problem::PoiseuilleFlow, x::Float64, y::Float64, timestep::Float64 = 0.0)
    return 1.0
end

function velocity(problem::PoiseuilleFlow, x::Float64, y::Float64, timestep::Float64 = 0.0)
    G = 1.0

    return [
        y * (problem.domain_size[2] - y) * (G / 2)
        0.0
    ]
end
function force(problem::PoiseuilleFlow, x::Float64, y::Float64, time::Float64 = 0.0)
    G = 1.0
    ν = viscosity(problem)

    return [
        ν * G
        0.0
    ]
end

function delta_x(problem::PoiseuilleFlow)
    return problem.domain_size[2] / problem.NY
end

has_external_force(problem::PoiseuilleFlow) = true

function apply_boundary_conditions_before!(q::Quadrature, problem::PoiseuilleFlow; time = t * Δt, f_new, f_old)
    for k = 1:size(f_old, 3)
        f_old[1, :, k] = f_old[problem.NX - 1, :, k]
        f_old[problem.NX, :, k]= f_old[2, :, k]
    end
end

"""
Apply a bottom and top wall bounce back boundary condition

The nodes at y = 1, NY are facing a stationary wall
The bounce back boundary condition will reflect any incoming particles
"""
function apply_boundary_conditions_after!(q::Quadrature, problem::PoiseuilleFlow; time = t * Δt, f_new, f_old)
    # Bottom wall
    y_idx = 1
    for f_idx = 1:size(f_new, 3)
        if q.abscissae[2, f_idx] > 0
            f_new[:, y_idx, f_idx] = f_old[:, y_idx, opposite(q, f_idx)]
        end
    end

    # Top wall
    y_idx = problem.NY
    for f_idx = 1:size(f_new, 3)
        if q.abscissae[2, f_idx] < 0
            f_new[:, y_idx, f_idx] = f_old[:, y_idx, opposite(q, f_idx)]
        end
    end
end
