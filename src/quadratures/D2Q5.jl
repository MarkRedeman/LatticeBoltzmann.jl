# https://github.com/grasingerm/lbxflow/blob/master/inc/lattice.jl
# http://mafija.fmf.uni-lj.si/seminar/files/2012_2013/Lattice_Boltzmann_method.pdf

# const original_order = [1, 2 ,3, 4]

# From: E3,1,5
struct D2Q5 <: Quadrature
    abscissae::Array{Int64, 2}
    weights::Array{Float64, 1}
    speed_of_sound_squared::Float64
end

D2Q5() = D2Q5(
    [
        0   1  -1   0   0
        0   0   0   1  -1
    ],
    [
        4/6, 1/12, 1/12, 1/12, 1/12
    ],
    6.0
)

order(q::D2Q5) = 3
Base.show(io::IO, q::D2Q5)= show(io, "D2Q5")
