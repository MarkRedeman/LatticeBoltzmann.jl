density(q::Quadrature, f::Array{Float64,N}) where {N} = sum(f, dims = N)
density(q::Quadrature, f::Array{Float64,1}) = sum(f)

function velocity!(q::Quadrature, f::Array{Float64,1}, ρ::Float64, u::Array{Float64,1})
    @inbounds for d = 1:dimension(q)
        u[d] = 0.0
        for idx = 1:length(f)
            u[d] += f[idx] * q.abscissae[d, idx]
        end
        u[d] / ρ
    end
    return
end

function pressure(
    q::Quadrature,
    f::Array{Float64},
    ρ::Float64,
    u::Array{Float64,1},
)::Float64
    a_2 = sum(f[idx] * hermite(Val{2}, q.abscissae[:, idx], q) for idx = 1:length(q.weights))
    D = dimension(q)

    p = (tr(a_2) - ρ * (u[1]^2 + u[2]^2 - D)) / D
    return p

    a_eq_2 = equilibrium_coefficient(Val{2}, q, ρ, u, 1.0)
    P = tr(a_eq_2) - ρ * (u[1]^2 + u[2]^2)

    D = dimension(q)
    E = 0.0
    @inbounds for idx = 1:length(f)
        E += f[idx] * (q.abscissae[1, idx]^2 + q.abscissae[2, idx]^2)
    end

    @show p q.speed_of_sound_squared * (E - ρ * (u[1]^2 + u[2]^2)) / D
    p = q.speed_of_sound_squared * (E - ρ * (u[1]^2 + u[2]^2)) / D

    return p
end
function momentum_flux(q::Quadrature, f::Array{Float64}, ρ::Float64, u::Array{Float64,1})
    D = dimension(q)
    P = zeros(D, D)
    @inbounds for x_idx = 1:D, y_idx = 1:D
        P[x_idx, y_idx] = sum(
            # Pressure tensor
            f[f_idx] *
            (q.abscissae[x_idx, f_idx] - u[x_idx]) *
            (q.abscissae[y_idx, f_idx] - u[y_idx])

            # Stress tensor
            # f[f_idx] * (q.abscissae[x_idx, f_idx]) * (q.abscissae[y_idx, f_idx])
            for f_idx = 1:length(f)
        ) #- ρ * u[x_idx] * u[y_idx]
    end

    return q.speed_of_sound_squared * P #- I * pressure(q, f, ρ, u)

    E = 0.0
    @inbounds for idx = 1:length(f)
        E += f[idx] * (q.abscissae[1, idx]^2 + q.abscissae[2, idx]^2)
    end

    p = q.speed_of_sound_squared * (E - ρ * (u[1]^2 + u[2]^2)) / D
    return p
end

# total energy density
function total_energy(q::Quadrature, f)
    return sum(
        [
            f[x, y, f_idx] * (q.abscissae[:, f_idx]' * q.abscissae[:, f_idx])
            for x = 1:size(f, 1), y = 1:size(f, 2), f_idx = 1:size(f, 3)
        ],
        dims = 3,
    )
end

# internal energy density
function internal_energy(q::Quadrature, f, ρ, u)
    return total_energy(q, f) - kinetic_energy(q, f, ρ, u)
end

# kinetic energy density
function kinetic_energy(q::Quadrature, f, ρ, u)
    return [
        ρ[x, y] * (u[x, y, 1] .^ 2 + u[x, y, 2] .^ 2)
        for x = 1:size(f, 1), y = 1:size(f, 2)
    ]
end

dimension(q::Q) where { Q <: Quadrature } = 2

function temperature(q::Quadrature, f, ρ, u)
    return pressure(q, f, ρ[:, :, 1], u) ./ ρ
    # return internal_energy(q, f, ρ, u) * (2 / dimension(q))
end

function temperature(q::Quadrature, f::Vector{Float64}, ρ::Float64, u::Vector{Float64})
    return pressure(q, f, ρ, u) ./ ρ
end

"""
Computes the deviatoric tensor σ

τ is the relaxation time such that ν = cs^2 τ
"""
function deviatoric_tensor(q::Quadrature, τ, f::Vector{Float64}, ρ::Float64, u::Vector{Float64})
    D = dimension(q)

    a_bar_2 = sum([f[idx] * hermite(Val{2}, q.abscissae[:, idx], q) for idx = 1:length(q.weights)])
    a_eq_2 = equilibrium_coefficient(Val{2}, q, ρ, u, 1.0)
    σ = (a_bar_2 - a_eq_2) / (1 + 1 / (2 * τ))
    return σ - I * tr(σ) / D
end
