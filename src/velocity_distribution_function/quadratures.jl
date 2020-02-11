# This file contains equilibria for specific quadratures

function _equilibrium(q::D2Q17, ρ, weight, u_dot_xi, u_squared, T, xi_squared)
    # Truncated upto order 2
    cs = q.speed_of_sound_squared
    D = dimension(q)
    # H_2_temperature = cs * (T .- 1) .* (cs * xi_squared - D)
    # H_3_temperature = 3.0 * cs * (T .- 1) * (cs * xi_squared - 2 - D)
    H_2_temperature = 0.0
    H_3_temperature = 0.0

    a_H_0 = 1.0
    a_H_1 = cs * u_dot_xi
    a_H_2 = cs^2 * (u_dot_xi .* u_dot_xi) .+ H_2_temperature .+ -cs * u_squared
    a_H_3 =
        cs * u_dot_xi .*
        (cs^2 * (u_dot_xi .* u_dot_xi) - 3 * cs * u_squared .+ H_3_temperature)
    return ρ .* weight .* (a_H_0 .+ a_H_1 .+ (1 / 2) * a_H_2 .+ (1 / 6) * a_H_3)
end

# function _equilibrium(q::D2Q17, ρ, weight, u_dot_xi, u_squared, T, xi_squared)
#     # Truncated upto order 2
#     cs = q.speed_of_sound_squared
#     D = dimension(q)
#     # H_2_temperature = (cs * T .- 1) .* (cs * xi_squared - D)
#     # H_3_temperature = 3.0 * (cs * T .- 1) * (cs * xi_squared - 2 - D)
#     H_2_temperature = 0.0
#     H_3_temperature = 0.0

#     a_H_0 = 1.0
#     a_H_1 = cs * u_dot_xi

#     a_H_2 = cs^2 * (u_dot_xi .* u_dot_xi) .+ H_2_temperature .+ - cs * u_squared
#     a_H_3 = cs * u_dot_xi .* (
#         cs^2 * (u_dot_xi .* u_dot_xi) -  3 * cs * u_squared .+ H_3_temperature
#     )
#     return ρ .* weight .* (
#         a_H_0 .+
#         a_H_1 .+
#         (1 / 2) * a_H_2 .+
#         (1 / 6) * a_H_3
#     )
# end

function _equilibrium(q::D2Q21, ρ, weight, u_dot_xi, u_squared, T, xi_squared)
    # Truncated upto order 2
    cs = q.speed_of_sound_squared
    D = dimension(q)
    # H_2_temperature = (cs * T .- 1) .* (cs * xi_squared - D)
    # H_3_temperature = 3.0 * (cs * T .- 1) * (cs * xi_squared - 2 - D)
    H_2_temperature = 0.0
    H_3_temperature = 0.0

    a_H_0 = 1.0
    a_H_1 = cs * u_dot_xi
    a_H_2 = cs^2 * (u_dot_xi .* u_dot_xi) .+ H_2_temperature .+ -cs * u_squared
    a_H_3 =
        cs * u_dot_xi .*
        (cs^2 * (u_dot_xi .* u_dot_xi) - 3 * cs * u_squared .+ H_3_temperature)
    return ρ .* weight .* (a_H_0 .+ a_H_1 .+ (1 / 2) * a_H_2 .+ (1 / 6) * a_H_3)
end

function equilibrium!(
    q::D2Q37,
    ρ::T,
    u::VT,
    temperature::T,
    f::VT,
) where {T <: Real, VT <: AbstractVector{T}}
    u_squared = zero(T)
    u_fourth = zero(T)
    for d in 1:dimension(q)
        u_squared += u[d] .^ 2
        u_fourth += u[d] .^ 4
    end

    for idx in 1:length(q.weights)
        u_dot_xi = q.abscissae[1, idx] .* u[1] .+ q.abscissae[2, idx] .* u[2]

        f[idx] = _equilibrium(
            q,
            ρ,
            q.weights[idx],
            u_dot_xi,
            u_squared,
            u_fourth,
            temperature,
            q.abscissae[1, idx]^2 + q.abscissae[2, idx]^2,
        )
    end

    return
end

function _equilibrium(
    q::D2Q37,
    ρ::T,
    weight::T,
    u_dot_xi::T,
    u_squared::T,
    u_fourth::T,
    temperature::T,
    xi_squared,
) where {T <: Real}
    cs = q.speed_of_sound_squared
    D = dimension(q)
    # H_2_temperature = (cs * temperature .- 1) .* (cs * xi_squared - D)
    # H_3_temperature = 3.0 * (cs * temperature .- 1) * (cs * xi_squared - 2 - D)
    H_2_temperature = zero(T)
    H_3_temperature = zero(T)

    a_H_0 = 1.0
    a_H_1 = cs * u_dot_xi
    a_H_2 = cs^2 * (u_dot_xi .* u_dot_xi) .+ H_2_temperature .+ -cs * u_squared
    a_H_3 =
        cs * u_dot_xi .*
        (cs^2 * (u_dot_xi .* u_dot_xi) - 3 * cs * u_squared .+ H_3_temperature)

    # H_4_temperature = 0.0
    a_H_4 = (cs^4 * u_dot_xi^4 - 6 * cs^3 * u_squared * u_dot_xi^2 + 3 * cs^2 * u_squared^2)

    return ρ .* weight .*
           (a_H_0 .+ a_H_1 .+ (1 / 2) * a_H_2 .+ (1 / 6) * a_H_3 .+ (1 / 24) * a_H_4)
end

# The D2Q4 lattice is unable to include temperature
pressure(q::D2Q4, f::VT, ρ::T, u::VT) where {T <: Real, VT <: AbstractVector{T}} = 1.0

function _equilibrium(
    q::D2Q4,
    ρ::T,
    weight::T,
    u_dot_xi::T,
    u_squared::T,
    temperature::T,
    xi_squared,
) where {T <: Real}
    cs = q.speed_of_sound_squared
    a_H_0 = one(T)
    a_H_1 = cs * u_dot_xi

    return ρ .* weight .* (a_H_0 .+ a_H_1)
end

function _equilibrium(
    q::D2Q5,
    ρ::T,
    weight::T,
    u_dot_xi::T,
    u_squared::T,
    temperature::T,
    xi_squared,
) where {T <: Real}
    cs = q.speed_of_sound_squared
    a_H_0 = one(T)
    a_H_1 = cs * u_dot_xi

    return ρ .* weight .* (a_H_0 .+ a_H_1)
end
