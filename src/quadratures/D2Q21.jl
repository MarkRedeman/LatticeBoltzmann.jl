# D2Q21
# See lectures in LBM

# Same as D2Q17, but vertical and horizontal vectors are extended by one

struct D2Q21{
    Abscissaes <: AbstractMatrix{Int64},
    Weights <: AbstractVector{Float64}
} <: Quadrature
    abscissae::Abscissaes
    weights::Weights
    speed_of_sound_squared::Float64
end
function D2Q21()
    cs = sqrt(2 / 3)
    w_0 = 91 / 324
    w_1 = 1 / 12
    w_2 = 2 / 27
    w_3 = 7 / 360
    w_4 = 1 / 432
    w_5 = 0 / 36
    w_6 = 1 / 1620

    # https://strathprints.strath.ac.uk/29195/1/Zhang_YH_Pure_Accuracy_analysis_of_high_order_lattice_Boltzmann_models_for_rarefied_gas_flows_21_Jan_2011.pdf
    w_0 = 91 / 324
    w_1 = 1 / 12
    w_2 = 2 / 27
    w_3 = 7 / 360
    w_4 = 1 / 432
    w_5 = 1 / 1620
    w_6 = 0
    cs = 3 / 2

    return D2Q21(
        [   #A #B          #C          #D          #E          #F
            0 1 -1 0 0 1 -1 1 -1 2 -2 0 0 2 -2 2 -2 3 -3 0 0 3 -3 -3 3
            0 0 0 1 -1 1 -1 -1 1 0 0 2 -2 2 -2 -2 2 0 0 3 -3 3 -3 3 -3
        ],
        [
            w_0,
            w_1,
            w_1,
            w_1,
            w_1,
            w_2,
            w_2,
            w_2,
            w_2,
            w_3,
            w_3,
            w_3,
            w_3,
            w_4,
            w_4,
            w_4,
            w_4,
            w_5,
            w_5,
            w_5,
            w_5,
            w_6,
            w_6,
            w_6,
            w_6,
        ],
        cs,
    )
end

order(q::D2Q21) = 7
Base.show(io::IO, q::D2Q21) = show(io, "D2Q21")
Base.string(q::D2Q21) = "D2Q21"


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
