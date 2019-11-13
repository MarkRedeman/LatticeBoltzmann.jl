using LinearAlgebra

abstract type Quadrature end

"""
Get the index of the abscissae pointing in the opposite direction of the given index

Here we assume that the abscissae of a quadrature are ordered such that each abscissae
follows its opposite abscissae
"""
function opposite(q::Quadrature, idx::Int64)
    if idx == 1
        return 1
    end
    if (mod(idx, 2) == 0)
        return idx + 1
    end
    return idx - 1
end

δ(α, β) = α == β ? 1 : 0
include("hermite-polynomials/hermite.jl")
include("hermite-polynomials/equilibrium-coefficients.jl")
include("hermite-polynomials/equilibrium.jl")
include("quadratures/quadrature.jl")
include("thermodynamics.jl")
