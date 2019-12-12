struct D2Q37 <: Quadrature
    abscissae::Array{Int64,2}
    weights::Array{Float64,1}
    speed_of_sound_squared::Float64
end

D2Q37() = D2Q37(
    [
        #  # group 2   # group 3   # group 4   # group 5               # group 6   # group 7   # group 8
        0 1 -1 0 0 1 -1 1 -1 2 -2 0 0 2 -2 -2 2 1 -1 1 -1 2 -2 2 -2 3 -3 0 0 3 -3 3 -3 1 -1 -1 1
        0 0 0 1 -1 1 -1 -1 1 0 0 2 -2 1 -1 1 -1 2 -2 -2 2 2 -2 -2 2 0 0 3 -3 1 -1 -1 1 3 -3 3 -3
    ],
    [
        0.23315066913235250228650

        # group 2
        0.10730609154221900241246
        0.10730609154221900241246
        0.10730609154221900241246
        0.10730609154221900241246

        # group 3
        0.05766785988879488203006
        0.05766785988879488203006
        0.05766785988879488203006
        0.05766785988879488203006

        # group 4
        0.01420821615845075026469
        0.01420821615845075026469
        0.01420821615845075026469
        0.01420821615845075026469

        # group 5
        0.00535304900051377523273
        0.00535304900051377523273
        0.00535304900051377523273
        0.00535304900051377523273
        0.00535304900051377523273
        0.00535304900051377523273
        0.00535304900051377523273
        0.00535304900051377523273

        # group 6
        0.00101193759267357547541
        0.00101193759267357547541
        0.00101193759267357547541
        0.00101193759267357547541

        # group 7
        0.00024530102775771734547
        0.00024530102775771734547
        0.00024530102775771734547
        0.00024530102775771734547

        # group 8
        0.00028341425299419821740
        0.00028341425299419821740
        0.00028341425299419821740
        0.00028341425299419821740
        0.00028341425299419821740
        0.00028341425299419821740
        0.00028341425299419821740
        0.00028341425299419821740
    ],
    # r =
    1.19697977039307435897239^2,
)

order(q::D2Q37) = 9
Base.show(io::IO, q::D2Q37) = show(io, "D2Q37")
function equilibrium!(q::D2Q37, ρ::Float64, u::Array{Float64,1}, T::Float64, f)
    u_squared = 0.0
    u_fourth = 0.0
    for d = 1:dimension(q)
        u_squared += u[d] .^ 2
        u_fourth += u[d] .^ 4
    end

    for idx = 1:length(q.weights)
        u_dot_xi = q.abscissae[1, idx] .* u[1] .+ q.abscissae[2, idx] .* u[2]

        f[idx] = _equilibrium(
            q,
            ρ,
            q.weights[idx],
            u_dot_xi,
            u_squared,
            u_fourth,
            T,
            q.abscissae[1, idx]^2 + q.abscissae[2, idx]^2,
        )
    end

    return
end

function _equilibrium(q::D2Q37, ρ, weight, u_dot_xi, u_squared, u_fourth, T, xi_squared)
    cs = q.speed_of_sound_squared
    D = dimension(q)
    # H_2_temperature = (cs * T .- 1) .* (cs * xi_squared - D)
    # H_3_temperature = 3.0 * (cs * T .- 1) * (cs * xi_squared - 2 - D)
    H_2_temperature = 0.0
    H_3_temperature = 0.0

    a_H_0 = 1.0
    a_H_1 = cs * u_dot_xi
    a_H_2 = cs^2 * (u_dot_xi .* u_dot_xi) .+ H_2_temperature .+ - cs * u_squared
    a_H_3 = cs * u_dot_xi .* (
        cs^2 * (u_dot_xi .* u_dot_xi) -  3 * cs * u_squared .+ H_3_temperature
    )

    # H_4_temperature = 0.0
    a_H_4 =  (
        cs^4 * u_dot_xi^4 - 6 * cs^3 * u_squared * u_dot_xi^2 + 3 * cs^2 * u_squared^2
    )

    return ρ .* weight .* (
        a_H_0 .+
        a_H_1 .+
        (1 / 2) * a_H_2 .+
        (1 / 6) * a_H_3 .+
        (1 / 24) * a_H_4
    )
end
