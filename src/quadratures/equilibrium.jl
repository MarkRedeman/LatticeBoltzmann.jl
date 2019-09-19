function equilibrium(
    q::Quadrature,
    ρ,#::Array{Float64, 2},
    u,#::Array{Float64, 3},
    T#::Array{Float64, 2}
)
    f = zeros(size(ρ,1), size(ρ,2), length(q.weights));

    u_squared = u[:, :, 1].^2 + u[:, :, 2].^2
    for idx = 1:length(q.weights)
        cs = q.abscissae[1, idx] .* u[:, :, 1] .+ q.abscissae[2, idx] .* u[:, :, 2]

        f[:, :, idx] = _equilibrium(q, ρ, q.weights[idx], cs, u_squared, T, q.abscissae[1, idx]^2 + q.abscissae[2, idx]^2)
    end

    return f
end
function equilibrium!(
    q::Quadrature,
    ρ::Array{Float64, 2},
    u::Array{Float64, 3},
    T::Array{Float64, 2},
    f
)
    for x = 1 : size(ρ, 1), y = 1 : size(ρ, 2)
        u_squared = u[x, y, 1].^2 + u[x, y, 2].^2
        for idx = 1:length(q.weights)
            cs = q.abscissae[1, idx] .* u[x, y, 1] .+ q.abscissae[2, idx] .* u[x, y, 2]

            f[x, y, idx] = _equilibrium(
                q,
                ρ[x, y],
                q.weights[idx],
                cs,
                u_squared,
                T[x, y],
                q.abscissae[1, idx]^2 + q.abscissae[2, idx]^2
            )
        end
    end
end

function equilibrium!(
    q::Quadrature,
    ρ::Array{Float64, 2},
    u::Array{Float64, 3},
    T::Float64,
    f
)
    for x = 1 : size(ρ, 1), y = 1 : size(ρ, 2)
        u_squared = u[x, y, 1].^2 + u[x, y, 2].^2
        for idx = 1:length(q.weights)
            cs = q.abscissae[1, idx] .* u[x, y, 1] .+ q.abscissae[2, idx] .* u[x, y, 2]

            f[x, y, idx] = _equilibrium(
                q,
                ρ[x, y],
                q.weights[idx],
                cs,
                u_squared,
                T,
                q.abscissae[1, idx]^2 + q.abscissae[2, idx]^2
            )
        end
    end
end

function equilibrium(
    q::Quadrature,
    ρ::Float64,
    u::Array{Float64, 1},
    T::Float64
)
    u_squared = sum(u.^2)
    f = zeros(length(q.weights));

    for idx = 1:length(q.weights)
        cs = q.abscissae[1, idx] .* u[1] .+ q.abscissae[2, idx] .* u[2]

        f[idx] = _equilibrium(
            q,
            ρ,
            q.weights[idx],
            cs,
            u_squared,
            T,
            q.abscissae[1, idx]^2 + q.abscissae[2, idx]^2
        )
    end
    return f
end


function equilibrium!(
    q::Quadrature,
    ρ::Float64,
    u::Array{Float64, 1},
    T::Float64,
    f
)
    u_squared = sum(u.^2)

    for idx = 1:length(q.weights)
        cs = q.abscissae[1, idx] .* u[1] .+ q.abscissae[2, idx] .* u[2]

        f[idx] = _equilibrium(
            q,
            ρ,
            q.weights[idx],
            cs,
            u_squared,
            T,
            q.abscissae[1, idx]^2 + q.abscissae[2, idx]^2
        )
    end
end

# Truncated upto order 2
function _equilibrium(q::Quadrature, ρ, weight, u_dot_xi, u_squared, T, xi_squared)
    cs = q.speed_of_sound_squared
    # a = 0.0
    a_H_0 = 1.0
    a_H_1 = cs * u_dot_xi

    H_2_temperature = (cs * T .- 1) .* (xi_squared * cs - dimension(q))
    # H_2_temperature = 0.0
    a_H_2 = cs^2 * (u_dot_xi .* u_dot_xi) .+ H_2_temperature .+ - cs * u_squared
    return ρ .* weight .* (
        a_H_0 .+
        a_H_1 .+
        (1 / 2) * a_H_2
    )
end

function _equilibrium(q::D2Q4, ρ, weight, u_dot_xi, u_squared, T, xi_squared)
    cs = q.speed_of_sound_squared
    a_H_0 = 1.0
    a_H_1 = cs * u_dot_xi

    return ρ .* weight .* (a_H_0 .+ a_H_1)
end

function _equilibrium(q::D2Q5, ρ, weight, u_dot_xi, u_squared, T, xi_squared)
    cs = q.speed_of_sound_squared
    a_H_0 = 1.0
    a_H_1 = cs * u_dot_xi

    return ρ .* weight .* (a_H_0 .+ a_H_1)
end

function _equilibrium(q::D2Q17, ρ, weight, u_dot_xi, u_squared, T, xi_squared)
    # Truncated upto order 2
    cs = q.speed_of_sound_squared
    D = dimension(q)
    H_2_temperature = (cs * T .- 1) .* (cs * xi_squared - D)
    H_3_temperature = 3.0 * (cs * T .- 1) * (cs * xi_squared -2 - D)
    H_3_temperature = 0.0

    a_H_0 = 1.0
    a_H_1 = cs * u_dot_xi
    a_H_2 = cs^2 * (u_dot_xi .* u_dot_xi) .+ H_2_temperature .+ - cs * u_squared
    a_H_3 = cs * u_dot_xi .* (
        cs^2 * (u_dot_xi .* u_dot_xi) -  3 * cs * u_squared .+ H_3_temperature
    )
    return ρ .* weight .* (
        a_H_0 .+
        a_H_1 .+
        (1 / 2) * a_H_2 .+
        (1 / 6) * a_H_3
    )
end


function hermite_equilibrium(q::Quadrature, f)
    ρ = density(q, f)
    j = momentum(q, f)
    u = j ./ ρ
    T = temperature(q, f, ρ, u)

    ξ = q.abscissae
    u_dot_ξ = dot(u, ξ)

    return ρ + ρ * u' * u + ρ * T * u' * u
end

function hermite_first_nonequilibrium(q::Quadrature, f)
end
