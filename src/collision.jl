using TimerOutputs

abstract type CollisionModel end
struct SRT <: CollisionModel
    τ::Float64
end

function collide(collision_model::SRT, q, f_in)::Array{Float64, 3}
    τ = collision_model.τ

    # Density
    ρ = density(q, f_in)

    # Momentum
    j = momentum(q, f_in) 

    # Temperature
    T = 1.0
    # T = temperature(q, f_in, ρ, j ./ ρ)

    feq = equilibrium(q, ρ, j ./ ρ, T);

    return (1 - 1 / τ) * f_in + (1 / τ) * feq;
    # Since this is the place where we will have computed all
    # hermite coefficients, call process! here? (or at least a callback)
    # cb(f_out, moments)
end

struct SRT_Force <: CollisionModel
    τ::Float64
    force::Array{Float64, 3}
end

function collide(collision_model::SRT_Force, q, f_in)::Array{Float64, 3}
    τ = collision_model.τ

    # Density
    ρ = density(q, f_in)

    # Momentum
    j = momentum(q, f_in)

    # Temperature
    T = 1.0
    # T = temperature(q, f_in, ρ, j ./ ρ)

    feq = equilibrium(
        q,
        ρ,
        j ./ ρ .+ τ * (0.0001) * collision_model.force,
        T
    );

    return (1 - 1 / τ) * f_in + (1 / τ) * feq;
    # Since this is the place where we will have computed all
    # hermite coefficients, call process! here? (or at least a callback)
    # cb(f_out, moments)
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

function collide(collision_model::SRT_Force, q, f_in, t)
    τ = collision_model.τ

    # Density
    ρ = density(q, f_in)

    # Momentum
    j = momentum(q, f_in)

    # Temperature
    T = 1.0
    # T = temperature(q, f_in, ρ, j ./ ρ)

    Δx = 1.0 / problem.N
    Δt = 1.0 / problem.N_iter
    force = (Δt^2 / Δx) * (Δt / Δx) * [
        TaylorGreen.force(problem.model, x, y) for x = linspace(0.0, 1.0, problem.N), y = linspace(0.0, 1.0, problem.N)
    ]

    feq = equilibrium(
        q,
        ρ,
        j ./ ρ + τ * force,
        T
    );

    f_out = (1 - 1 / τ) * f_in + (1 / τ) * feq;


    # Since this is the place where we will have computed all
    # hermite coefficients, call process! here? (or at least a callback)
    # cb(f_out, moments)

    return f_out
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
