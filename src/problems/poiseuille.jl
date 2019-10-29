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
    u_max = 0.1 / scale

    Re = NY * u_max / ν
    @show u_max, Re, ν

    # @show (ν - Δt / 2)^2
    return PoiseuilleFlow(
        1.0,
        u_max,
        ν,
        # NX,
        5,
        NY,
        1.0,
        domain_size
    )
end

function range(problem::PoiseuilleFlow)
    Δx = problem.domain_size[1] / problem.NX
    Δy = problem.domain_size[2] / problem.NY

    x_range = range(Δx / 2, problem.domain_size[1] - Δx / 2, length = problem.NX)
    y_range = range(Δy / 2, problem.domain_size[1] - Δy / 2, length = problem.NY)

    return x_range, y_range
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
    # Don't include the two boundary nodes
    return problem.domain_size[2] * (1 / (problem.NY) )
end

has_external_force(problem::PoiseuilleFlow) = true

# Temporary
function is_fluid(problem::PoiseuilleFlow, x::Int64, y::Int64)
    return true
    if y == 1
        return false
    end
    if y == problem.NY
        return false
    end
    return true
end

function apply_boundary_conditions_before!(q::Quadrature, problem::PoiseuilleFlow; time = t * Δt, f_new, f_old)
    for k = 1:size(f_old, 3)
        f_old[1, :, k] = f_old[problem.NX - 1, :, k]
        f_old[problem.NX, :, k]= f_old[2, :, k]
    end
end

function apply_boundary_conditions_after!(q::Quadrature, problem::PoiseuilleFlow; time = t * Δt, f_new, f_old)
    for x_idx = 1 : problem.NX
        # Bottom wall
        y_idx = 1
        # for f_idx = 1:size(f_new, 3)
        for f_idx = 1:size(f_new, 3)
            if q.abscissae[2, f_idx] > 0
                f_new[x_idx, y_idx, f_idx] = f_old[x_idx, y_idx, opposite(q, f_idx)]
            end
        end
        # for f_idx = [4, 5, 6]
        #     f_new[x_idx, y_idx, f_idx] = f_old[x_idx, y_idx, opposite(q, f_idx)]
        # end

        # Top wall
        y_idx = problem.NY

        for f_idx = 1:size(f_new, 3)
            if q.abscissae[2, f_idx] < 0
                f_new[x_idx, y_idx, f_idx] = f_old[x_idx, y_idx, opposite(q, f_idx)]
            end
        end
        # for f_idx = 1:size(f_new, 3)
        # for f_idx = [2, 8, 9]
        #     f_new[x_idx, y_idx, f_idx] = f_old[x_idx, y_idx, opposite(q, f_idx)]
        # end
    end
end
