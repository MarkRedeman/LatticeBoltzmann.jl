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

function range(problem::PoiseuilleFlow)
    x_range = range(0.0, problem.domain_size[1], length=problem.NX + 1)
    y_range = range(0.0, problem.domain_size[2], length=problem.NY + 1)

    return x_range, y_range
end
function PoiseuilleFlow(
    ν = 1.0 / 6.0, scale = 2, NX = 8 * scale + 0, NY = NX, domain_size = (2pi, 2pi); static = true
)
    u_max = 0.025 / scale
    Re = NX * u_max / ν
    @show u_max, Re
    return PoiseuilleFlow(
        1.0,
        u_max,
        ν,
        NX,
        NY + 2,
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
    ν = viscosity(problem)

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

apply_boundary_conditions!(q, problem; time = t * Δt, f_new, f_old) = apply_boundary_conditions!(q, problem, f_old, f_new, time = time)
function apply_boundary_conditions!(q::Quadrature, problem::PoiseuilleFlow, f_in, f_out; time = 0.0)
    for y_idx = [1, problem.NY], x_idx = 1 : problem.NX
        f = copy(f_in[x_idx, y_idx, :])
        for f_idx = 1:size(f_out, 3)
            f_in[x_idx, y_idx, f_idx] = f[opposite(q, f_idx)]
        end
    end
end
