export SRT, TRT, MRT
export CollisionModel

abstract type CollisionModel end
struct SRT <: CollisionModel
    τ::Float64
end

struct TRT <: CollisionModel
    τ_symmetric
    τ_asymmetric
end
# Magic parameter
function TRT(τ_symmetric, Δt, ν)
    # \Nabla = (τ_symmetric / Δt - .5) * (τ_asymmetric / Δt - .5)
    # should equal 1 / 4
end

struct SRT_Force{T} <: CollisionModel
    τ::Float64
    # force::Array{Float64, 3}
    force::T#::Array{Any, 2}
end
SRT(τ, force) = SRT_Force(τ, force)


function CollisionModel(
    cm::Type{<:CollisionModel},
    q::Quadrature,
    problem::InitialValueProblem
)
    τ = q.speed_of_sound_squared * lattice_viscosity(problem) + 0.5

    if has_external_force(problem)
        force_field = (x_idx, y_idx, t) -> lattice_force(problem, x_idx, y_idx, t)
        return SRT(τ, force_field)
    end

    return SRT(τ)
end

collide!(c, q::Quadrature; time, f_new, f_old, problem) = collide!(c, q, f_old, f_new, time = time, problem = problem)

function collide!(collision_model::SRT, q::Quadrature, f_in, f_out; time = 0.0, problem = nothing)
    τ = collision_model.τ

    feq = Array{Float64}(undef, size(f_in, 3))
    f = Array{Float64}(undef, size(f_in, 3))
    u = zeros(dimension(q))
   
    @inbounds for x = 1 : size(f_in, 1), y = 1 : size(f_in, 2)
        @inbounds for f_idx = 1 : size(f_in, 3)
            f[f_idx] = f_in[x, y, f_idx]
        end

        ρ = density(q, f)

        # Momentum
        velocity!(q, f, ρ, u)

        # Temperature
        # T = temperature(q, f, ρ, u)
        T = 1.0

        equilibrium!(q, ρ, u, T, feq);

        @inbounds for f_idx = 1 : size(f_in, 3)
            f_out[x, y, f_idx] = (1 - 1 / τ) * f[f_idx] + (1 / τ) * feq[f_idx]
        end
    end
    return
end

function collide!(collision_model::SRT_Force, q::Quadrature, f_in, f_out; time = 0.0, problem = nothing)
    τ = collision_model.τ

    feq = Array{Float64}(undef, size(f_in, 3))
    f = Array{Float64}(undef, size(f_in, 3))
    u = zeros(dimension(q))
    F = zeros(dimension(q))

    @inbounds for x = 1 : size(f_in, 1), y = 1 : size(f_in, 2)
        @inbounds for f_idx = 1 : size(f_in, 3)
            f[f_idx] = f_in[x, y, f_idx]
        end

        ρ = density(q, f)

        # Momentum
        velocity!(q, f, ρ, u)

        # Temperature
        # T = temperature(q, f, ρ, u)
        T = 1.0

        F .= collision_model.force(x, y, time)

        equilibrium!(q, ρ, u + τ * F, T, feq);


        @inbounds for f_idx = 1 : size(f_in, 3)
            f_out[x, y, f_idx] = (1 - 1 / τ) * f[f_idx] + (1 / τ) * feq[f_idx];
        end
    end
    return
end

# function collide(collision_model::TRT, quadrature, f_in)
#     τ_s = collision_model.τ_symmetric
#     τ_a = collision_model.τ_asymmetric

#     # TODO
#     return f_in
# end
