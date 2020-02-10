abstract type Quadrature end

dimension(q::Q) where { Q <: Quadrature } = 2

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
include("hermite_polynomials/hermite.jl")
include("hermite_polynomials/equilibrium_coefficients.jl")
include("hermite_polynomials/equilibrium.jl")
include("thermodynamics.jl")

# Specific velocity sets
include("quadratures/D2Q4.jl")
include("quadratures/D2Q5.jl")
include("quadratures/D2Q9.jl")
include("quadratures/D2Q13.jl")
include("quadratures/D2Q17.jl")
include("quadratures/D2Q21.jl")
include("quadratures/D2Q37.jl")

const Quadratures = (
    D2Q4 = D2Q4(),
    D2Q5 = D2Q5(),
    D2Q9 = D2Q9(),
    D2Q13 = D2Q13(),
    D2Q17 = D2Q17(),
    D2Q21 = D2Q21(),
    D2Q37 = D2Q37(),
)

