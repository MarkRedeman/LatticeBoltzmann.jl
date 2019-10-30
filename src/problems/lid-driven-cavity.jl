export LidDrivenCavityFlow

struct LidDrivenCavityFlow <: InitialValueProblem
    rho_0::Float64
    u_max::Float64
    ν::Float64
    NX::Int64
    NY::Int64
    domain_size::Tuple{Float64, Float64}
end

function LidDrivenCavityFlow(
    ν = 1.0 / 6.0, scale = 2, NX = 16 * scale, NY = NX, domain_size = (1.0, 1.0)
)
    # # Since our velocity profile is constant over x it is sufficient to only take
    # # 3 nodes in the x direction
    # NX = 5
    u_max = 0.01 / scale

    return LidDrivenCavityFlow(1.0, u_max, ν, NX, NY, domain_size)
end

function density(q::Quadrature, problem::LidDrivenCavityFlow, x::Float64, y::Float64, timestep::Float64 = 0.0)
    return 1.0
end

function pressure(q::Quadrature, problem::LidDrivenCavityFlow, x::Float64, y::Float64, timestep::Float64 = 0.0)
    return 1.0
end

function velocity(problem::LidDrivenCavityFlow, x::Float64, y::Float64, timestep::Float64 = 0.0)
    return [
        0.0
        0.0
    ]
end

"""
Apply a bottom and top wall bounce back boundary condition

The nodes at y = 1, NY are facing a stationary wall
The bounce back boundary condition will reflect any incoming particles
"""
function apply_boundary_conditions_after!(q::Quadrature, problem::LidDrivenCavityFlow; time = t * Δt, f_new, f_old)
    bounce_back_from_left!(q, f_new, f_old, 1, :)
    bounce_back_from_right!(q, f_new, f_old, problem.NX, :)
    bounce_back_from_bottom!(q, f_new, f_old, :, 1)

    bounce_back_from_top!(q, f_new, f_old, :, problem.NY, [problem.u_max, 0])
end
