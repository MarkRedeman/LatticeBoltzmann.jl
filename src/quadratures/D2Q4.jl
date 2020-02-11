# https://github.com/grasingerm/lbxflow/blob/master/inc/lattice.jl
# http://mafija.fmf.uni-lj.si/seminar/files/2012_2013/Lattice_Boltzmann_method.pdf

# const original_order = [1, 2 ,3, 4]

struct D2Q4{
    T <: Real,
    Abscissaes <: AbstractMatrix{<:Integer},
    Weights <: AbstractVector{<:Real},
} <: Quadrature
    abscissae::Abscissaes
    weights::Weights
    speed_of_sound_squared::T
end

D2Q4() = D2Q4(
    [
        1 0 -1 0
        0 1 0 -1
    ],
    [1 / 4, 1 / 4, 1 / 4, 1 / 4],
    2.0,
)

order(q::D2Q4) = 3
function opposite(q::D2Q4, idx::Int64)
    if idx <= 2
        return idx + 2
    end
    return idx - 2
end

Base.show(io::IO, q::D2Q4) = show(io, "D2Q4")
Base.string(q::D2Q4) = "D2Q4"
