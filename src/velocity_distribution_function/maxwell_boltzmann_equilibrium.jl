function equilibrium(
    q::Quadrature,
    ρ::T,
    u::VT,
    temperature::T,
) where {T <: Real, VT <: AbstractVector{T}}
    f = zeros(length(q.weights))
    hermite_based_equilibrium!(q, ρ, u, temperature, f)
    return f
end

function equilibrium!(
    q::Quadrature,
    ρ::T,
    u::VT,
    temperature::T,
    f,
) where {T <: Real, VT <: AbstractVector{T}}
    u_squared = 0.0
    for d in 1:dimension(q)
        u_squared += u[d] .^ 2
    end

    for idx in 1:length(q.weights)
        u_dot_xi = q.abscissae[1, idx] .* u[1] .+ q.abscissae[2, idx] .* u[2]

        f[idx] = _equilibrium(
            q,
            ρ,
            q.weights[idx],
            u_dot_xi,
            u_squared,
            temperature,
            q.abscissae[1, idx]^2 + q.abscissae[2, idx]^2,
        )
    end

    return
end

# Truncated upto order 2
# function _equilibrium(q::Quadrature, ρ::T, weight::T, u_dot_xi::T, u_squared::T, temperature::T, xi_squared::T) where { T <: Real }
function _equilibrium(
    q::Quadrature,
    ρ,
    weight,
    u_dot_xi,
    u_squared,
    temperature,
    xi_squared,
)
    cs = q.speed_of_sound_squared
    # a = 0.0
    a_H_0 = 1.0
    a_H_1 = cs * u_dot_xi

    # H_2_temperature = cs * ( temperature .- 1) .* (xi_squared * cs - dimension(q))
    H_2_temperature = 0.0 #zero(T)

    a_H_2 = cs^2 * (u_dot_xi .* u_dot_xi) .+ H_2_temperature .+ -cs * u_squared
    # return ρ .* weight .+ weight * (
    #     a_H_1 .+
    #     (1 / 2) * a_H_2
    # )
    return ρ .* weight .* (a_H_0 .+ a_H_1 .+ (1 / 2) * a_H_2)
end
