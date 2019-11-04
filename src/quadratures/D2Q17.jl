# D2Q17
# See lectures in LBM

# Picture is D2Q9, but with each vector extended to a 3 x 3 box

# Contains description:
# Multiple-relaxation time model for the correct thermodynamic equations
#
# https://www.researchgate.net/publication/23315587_Multiple-relaxation_time_model_for_the_correct_thermodynamic_equations
#
# https://www.researchgate.net/profile/Zhaoli_Guo2/publication/23315587_Multiple-relaxation_time_model_for_the_correct_thermodynamic_equations/links/0046352b40f45699f8000000.pdf

function generate_quadratures()
    c = 1.
    rest = [0, 0]
    group_1 = [[cos(i - 1)π  / 2, sin(i - 1)π / 2] * c for i = 1 : 4]
    group_2 = [[cos(2i - 9)π  / 4, sin(2i - 9)π / 4] * sqrt(2) * c for i = 5 : 8]
    group_3 = [[cos(i - 1)π  / 2, sin(i - 1)π / 2] * 2c for i = 9 : 12]
    group_4 = [[cos(2i - 9)π  / 4, sin(2i - 9)π / 4] * sqrt(2) * c for i = 13 : 16]

    return [
        rest,
        group_1...,
        group_2...,
        group_3...,
        group_4...,
    ]
end


struct D2Q17 <: Quadrature
    abscissae::Array{Int64, 2}
    weights::Array{Float64, 1}
    speed_of_sound_squared::Float64

end
function D2Q17()
    sq = sqrt(193)
    w_0 = (575  + 193sq) / 8100
    w_1 = (3355 - 91sq)  / 18000
    w_2 = (655  + 17sq)  / 27000
    w_3 = (685  - 49sq)  / 54000
    w_4 = (1445 - 101sq) / 162000

    return D2Q17(
        [   #A #B          #C          #D          #D
            0  1 -1  0  0  1 -1  1 -1  2 -2  2 -2  3 -3  0  0
            0  0  0  1 -1  1 -1 -1  1  2 -2 -2  2  0  0  3 -3
            ],
        [
            w_0, w_1, w_1, w_1, w_1, w_2, w_2, w_2, w_2, w_3, w_3, w_3, w_3, w_4, w_4, w_4, w_4
        ],
        (125 + 5 * sqrt(193)) / 72
    )
end

order(q::D2Q17) = 7
Base.show(io::IO, q::D2Q17)= show(io, "D2Q17")
