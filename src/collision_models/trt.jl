struct TRT{Force, T <: Real} <: CollisionModel
    τ_symmetric::T
    τ_asymmetric::T
    force::Force
end
TRT(τ_a, τ_s) = TRT(τ_s, τ_a, nothing)

function CollisionModel(
    cm::Type{<:TRT},
    q::Quadrature,
    problem::FluidFlowProblem;
    Λ = 1 / 4,
)
    τ = q.speed_of_sound_squared * lattice_viscosity(problem) + 0.5
    τ_ = 0.5 + Λ / (τ - 0.5)

    force_field = has_external_force(problem) ?
        (x_idx, y_idx, t) -> lattice_force(problem, x_idx, y_idx, t) : nothing

    return TRT(τ, τ_, force_field)
end

"""
Acts like a factory for the TRT operator where the magic parameter Λ can be set

According to [TODO cite] Λ has the following effect:
- Λ = 1 / 12 # Cancels third-order spatial error, leading to optimal results for
  pure advection problems
- Λ = 1 / 6 # Cancels fourth-order spatial error, providing the most accurate
  results for the pure diffusion equation.
- Λ = 3 / 16 # Boundary wall location implemented via bounce-back for the
  Poiseuille flow exactly in the middle between horizontal walls and fluid nodes.
- Λ = 1 / 4 # Provides the most stable simulations
"""
struct TRT_Λ{T <: Real}
    Λ::T
end
function CollisionModel(cm::TRT_Λ, q::Quadrature, problem::FluidFlowProblem)
    CollisionModel(TRT, q, problem, Λ = cm.Λ)
end

function collide!(
    collision_model::TRT{Force},
    q::Quadrature,
    f_in::Populations,
    f_out::Populations;
    time = 0.0,
) where {Force, T <: Real, Populations <: AbstractArray{T, 3}}
    τ_s = collision_model.τ_symmetric
    τ_a = collision_model.τ_asymmetric

    feq = Array{T}(undef, size(f_in, 3))
    f = Array{T}(undef, size(f_in, 3))
    u = zeros(dimension(q))
    if !(Force <: Nothing)
        F = zeros(dimension(q))
    end

    nx, ny, nf = size(f_in)
    @inbounds for x = 1:nx, y = 1:ny
        @inbounds for f_idx = 1:nf
            f[f_idx] = f_in[x, y, f_idx]
        end

        ρ = density(q, f)

        # Momentum
        velocity!(q, f, ρ, u)

        # Temperature
        # temperature = temperature(q, f, ρ, u)
        temperature = 1.0

        if Force <: Nothing
            equilibrium!(q, ρ, u, temperature, feq)
        else
            F .= collision_model.force(x, y, time)

            equilibrium!(q, ρ, u + τ_a * F, temperature, feq)
        end

        @inbounds for f_idx = 1:nf
            opposite_idx = opposite(q, f_idx)
            feq_symmetric = 0.5 * (feq[f_idx] + feq[opposite_idx])
            feq_asymmetric = 0.5 * (feq[f_idx] - feq[opposite_idx])
            f_symmetric = 0.5 * (f[f_idx] + f[opposite_idx])
            f_asymmetric = 0.5 * (f[f_idx] - f[opposite_idx])


            f_out[x, y, f_idx] =
                f[f_idx] + (
                    -(1 / τ_s) * (f_symmetric - feq_symmetric) -
                    (1 / τ_a) * (f_asymmetric - feq_asymmetric)
                )
        end
    end
    return
end
