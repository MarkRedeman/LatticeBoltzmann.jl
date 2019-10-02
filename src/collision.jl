using TimerOutputs

abstract type CollisionModel end
struct SRT <: CollisionModel
    τ::Float64
end

struct SRT_Force{T} <: CollisionModel
    τ::Float64
    # force::Array{Float64, 3}
    force::T#::Array{Any, 2}
end

using TimerOutputs
function collide!(collision_model::SRT, q::Quadrature, f_in, f_out; time = 0.0, problem = nothing)
    τ = collision_model.τ

    feq = Array{Float64}(undef, size(f_in, 3))
    f = Array{Float64}(undef, size(f_in, 3))
    u = zeros(dimension(q))
   
    @inbounds for x = 1 : size(f_in, 1), y = 1 : size(f_in, 2)
        if ! is_fluid(problem, x, y)
            continue
        end


        @inbounds for f_idx = 1 : size(f_in, 3)
            f[f_idx] = f_in[x, y, f_idx]
        end

        ρ = density(q, f)

        # Momentum
        lbm.velocity!(q, f, ρ, u)

        # Temperature
        T = temperature(q, f, ρ, u)
        # T = 1.0 / q.speed_of_sound_squared

        equilibrium!(q, ρ, u, T, feq);


        @inbounds for f_idx = 1 : size(f_in, 3)
            f_out[x, y, f_idx] = (1 - 1 / τ) * f[f_idx] + (1 / τ) * feq[f_idx];
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
        if ! is_fluid(problem, x, y)
            continue
        end

        @inbounds for f_idx = 1 : size(f_in, 3)
            f[f_idx] = f_in[x, y, f_idx]
        end

        ρ = density(q, f)

        # Momentum
        lbm.velocity!(q, f, ρ, u)

        # Temperature
        T = temperature(q, f, ρ, u)
        # T = 1.0 / q.speed_of_sound_squared

        F .= τ * collision_model.force(x, y, time)

        equilibrium!(q, ρ, u + F, T, feq);


        @inbounds for f_idx = 1 : size(f_in, 3)
            f_out[x, y, f_idx] = (1 - 1 / τ) * f[f_idx] + (1 / τ) * feq[f_idx];
        end
    end
    return
end

function collide(collision_model::SRT, q, f_in; time = 0.0)::Array{Float64, 3}

    τ = collision_model.τ

    feq = Array{Float64}(undef, size(f_in, 3))
    f = Array{Float64}(undef, size(f_in, 3))
    j = zeros(dimension(q))
    f_out = copy(f_in)
    @inbounds for x = 1 : size(f_in, 1), y = 1 : size(f_in, 2)
        @inbounds for f_idx = 1 : size(f_in, 3)
            f[f_idx] = f_in[x, y, f_idx]
        end

        # Density
        ρ = density(q, f)

        # Momentum
        lbm.momentum!(q, f, j)

        # Temperature
        T = temperature(q, f, ρ, j ./ ρ)

        equilibrium!(q, ρ, j / ρ, T, feq);


        f_out[x, y, :] .= (1 - 1 / τ) * f + (1 / τ) * feq;
        # @inbounds for f_idx = 1 : size(f_in, 3)
        #     f_out[x, y, f_idx] = (1 - 1 / τ) * f[x, y, f_idx] + (1 / τ) * feq[f_idx];
        # end
    end
    return f_out
    τ = collision_model.τ

    # Density
    ρ = density(q, f_in)

    # Momentum
    j = momentum(q, f_in) 

    # Temperature
    # T = 1.0
    T = temperature(q, f_in, ρ, j ./ ρ)

    return srt_collision(
        q,
        f_in,
        τ,
        ρ,
        j ./ ρ,
        T
    )
end

function collide(collision_model::SRT_Force, q, f_in, time = 0.0)::Array{Float64, 3}
    τ = collision_model.τ

    # Density
    ρ = density(q, f_in)

    # Momentum
    j = momentum(q, f_in)

    # Temperature
    # T = 1.0
    T = temperature(q, f_in, ρ, j ./ ρ)

    return srt_collision(
        q,
        f_in,
        τ,
        ρ,
        j ./ ρ .+ τ * collision_model.force(time),
        T
    )
end

function srt_collision(q, f_in, τ, ρ, u, T)
    feq = equilibrium(q, ρ, u, T);

    return (1 - 1 / τ) * f_in + (1 / τ) * feq;
end


function collide_3(collision_model::SRT, q, f_in, f_out, feq)
    feq = Array{Float64, 3}(undef, size(f_in,1), size(f_in,2), length(q.weights));
    τ = collision_model.τ

    # Density
    ρ = [density(q, f_in[x, y, :]) for x = 1 : size(f_in, 1), y = 1 : size(f_in, 2)]

    # Momentum
    j = momentum(q, f_in)

    # Temperature
    T = 1.0
    # T = temperature(q, f_in, ρ, j ./ ρ)

    equilibrium!(q, ρ, j ./ ρ, T, feq);

    f_out = (1 - 1 / τ) * f_in + (1 / τ) * feq;

    return
    # @show to

    # Since this is the place where we will have computed all
    # hermite coefficients, call process! here? (or at least a callback)
    # cb(f_out, moments)
end

function collide_2(collision_model::SRT, q, f_in, f_out, to = TimerOutput())
    τ = collision_model.τ

    feq = Array{Float64}(undef, size(f_in, 3))
    f = Array{Float64}(undef, size(f_in, 3))
    j = zeros(size(f_in, 3))
    @inbounds for x = 1 : size(f_in, 1), y = 1 : size(f_in, 2)
        f = f_in[x, y, :]

        # Density
        ρ = density(q, f)

        # Momentum
        lbm.momentum!(q, f, j)

        # Temperature
        T = 1.0
        # T = temperature(q, f_in, ρ, j ./ ρ)

        equilibrium!(q, ρ, j / ρ, T, feq);

        f_out[x, y, :] = (1 - 1 / τ) * f + (1 / τ) * feq;
        # for f_idx = 1 : size(f_in, 3)
        #     f_out[x, y, f_idx] = (1 - 1 / τ) * f[f_idx] + (1 / τ) * feq[f_idx];
        # end
    end
    return
end

struct TRT <: CollisionModel
    τ_symmetric
    τ_asymmetric
end

function collide(collision_model::TRT, quadrature, f_in)
    τ_s = collision_model.τ_symmetric
    τ_a = collision_model.τ_asymmetric

    # TODO
    return f_in
end
