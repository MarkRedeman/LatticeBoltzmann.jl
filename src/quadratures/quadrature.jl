using LinearAlgebra

abstract type Quadrature end
# abstract type Quadrature{Thermal} end
abstract type Lattice end

DdQq(d, q) = "HOI"

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


include("hermite.jl")
include("D2Q4.jl")
include("D2Q5.jl")
include("D2Q9.jl")
include("D2Q17.jl")
include("thermodynamics.jl")
include("equilibrium.jl")
