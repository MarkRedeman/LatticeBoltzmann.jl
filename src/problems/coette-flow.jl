export CouetteFlow

import Base: range

struct CouetteFlow <: lbm.InitialValueProblem
    rho_0::Float64
    u_max::Float64
    ν::Float64
    NX::Int64
    NY::Int64
    k::Float64
    domain_size::Tuple{Float64, Float64}
end

function CouetteFlow(
    ν = 1.0 / 6.0, scale = 2, NX = 8 * scale + 0, NY = NX, domain_size = (1.0, 1.0); static = true
)
    u_max = 0.0125 / scale
    Re = NX * u_max / ν
    @show u_max, Re
    Δt = delta_t(CouetteFlow(1.0, u_max, ν, NX, NY + 2, 1.0, domain_size))
    @show (ν - Δt / 2)^2
    return CouetteFlow(
        1.0,
        u_max,
        ν,
        NX,
        NY + 2,
        1.0,
        domain_size
    )
end

function range(problem::CouetteFlow)
    x_range = range(0.0, problem.domain_size[1], length=problem.NX + 1)

    # do not add 1 since the boundary is inclusive
    y_range = range(0.0, problem.domain_size[2], length=problem.NY)


    # Accomodate for the wall at the midway boundary
    # y_range = range(0.0, problem.domain_size[2], length=problem.NY) # did not add 1
    # y_range = range(y_range[2], y_range[problem.NY], length=problem.NY - 2)

    return x_range, y_range
end

function density(q::Quadrature, problem::CouetteFlow, x::Float64, y::Float64, timestep::Float64 = 0.0)
    return 1.0
end

function pressure(q::Quadrature, problem::CouetteFlow, x::Float64, y::Float64, timestep::Float64 = 0.0)
    return 1.0
end

function velocity(problem::CouetteFlow, x::Float64, y::Float64, timestep::Float64 = 0.0)
    G = 1.0

    return [
        y * (problem.domain_size[2] - y) * (G / 2)
        0.0
    ]
end
function force(problem::CouetteFlow, x::Float64, y::Float64, time::Float64 = 0.0)
    G = 1.0
    ν = viscosity(problem)

    return [
        ν * G
        0.0
    ]
end

function delta_x(problem::CouetteFlow)
    # Don't include the two boundary nodes
    return problem.domain_size[2] * (1 / (problem.NY - 2) )
end

has_external_force(problem::CouetteFlow) = true

# Temporary
function is_fluid(problem::CouetteFlow, x::Int64, y::Int64)
    if y == 1
        return false
    end
    if y == problem.NY
        return false
    end
    return true
end

apply_boundary_conditions!(q, problem; time = t * Δt, f_new, f_old) = apply_boundary_conditions!(q, problem, f_old, f_new, time = time)
function apply_boundary_conditions!(q::Quadrature, problem::CouetteFlow, f_in, f_out; time = 0.0)

    f = Array{Float64}(undef, size(f_in, 3))

    for x_idx = 1 : problem.NX
        # Top wall
        y_idx = 1
        copy!(f, f_in[x_idx, y_idx, :])
        for f_idx = 1:size(f_out, 3)
            f_in[x_idx, y_idx, f_idx] = f[opposite(q, f_idx)]
        end

        # Bottom wall
        y_idx = problem.NY
        copy!(f, f_in[x_idx, y_idx, :])
        for f_idx = 1:size(f_out, 3)
            f_in[x_idx, y_idx, f_idx] = f[opposite(q, f_idx)]
        end
    end
end
