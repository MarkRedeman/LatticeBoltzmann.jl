using LinearAlgebra

abstract type Quadrature end
# abstract type Quadrature{Thermal} end
abstract type Lattice end

DdQq(d, q) = "HOI"

opposite(q::Quadrature, idx::Int64) = idx

include("hermite.jl")
include("D2Q4.jl")
include("D2Q5.jl")
include("D2Q9.jl")
include("D2Q17.jl")
include("thermodynamics.jl")
include("equilibrium.jl")
