# https://github.com/grasingerm/lbxflow/blob/master/inc/lattice.jl
# http://mafija.fmf.uni-lj.si/seminar/files/2012_2013/Lattice_Boltzmann_method.pdf

# const original_order = [1, 2 ,3, 4]

struct D2Q4{
    Abscissaes <: AbstractMatrix{Int64},
    Weights <: AbstractVector{Float64}
} <: Quadrature
    abscissae::Abscissaes
    weights::Weights
    speed_of_sound_squared::Float64
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

# The D2Q4 lattice is unable to include temperature
pressure(q::D2Q4, f::VT, ρ::Float64, u::VT) where { VT <: AbstractVector{Float64}} = 1.0

function _equilibrium(q::D2Q4, ρ, weight, u_dot_xi, u_squared, T, xi_squared)
    cs = q.speed_of_sound_squared
    a_H_0 = 1.0
    a_H_1 = cs * u_dot_xi

    return ρ .* weight .* (a_H_0 .+ a_H_1)
end
