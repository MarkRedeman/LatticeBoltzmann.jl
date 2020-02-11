# D2Q17
# See lectures in LBM

# Picture is D2Q9, but with each vector extended to a 3 x 3 box

# Contains description:
# Multiple-relaxation time model for the correct thermodynamic equations
#
# https://www.researchgate.net/publication/23315587_Multiple-relaxation_time_model_for_the_correct_thermodynamic_equations
#
# https://www.researchgate.net/profile/Zhaoli_Guo2/publication/23315587_Multiple-relaxation_time_model_for_the_correct_thermodynamic_equations/links/0046352b40f45699f8000000.pdf
struct D2Q17{
    T <: Real,
    Abscissaes <: AbstractMatrix{<:Integer},
    Weights <: AbstractVector{<:Real},
} <: Quadrature
    abscissae::Abscissaes
    weights::Weights
    speed_of_sound_squared::T
end
function D2Q17()
    sq = sqrt(193)
    w_0 = (575 + 193sq) / 8100
    w_1 = (3355 - 91sq) / 18000
    w_2 = (655 + 17sq) / 27000
    w_3 = (685 - 49sq) / 54000
    w_4 = (1445 - 101sq) / 162000

    return D2Q17(
        [   #A #B          #C          #D          #D
            0 1 -1 0 0 1 -1 1 -1 2 -2 2 -2 3 -3 0 0
            0 0 0 1 -1 1 -1 -1 1 2 -2 -2 2 0 0 3 -3
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
        ],
        (125 + 5 * sqrt(193)) / 72,
    )
end

order(q::D2Q17) = 7
Base.show(io::IO, q::D2Q17) = show(io, "D2Q17")
Base.string(q::D2Q17) = "D2Q17"
