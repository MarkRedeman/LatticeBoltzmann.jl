# Specific quadrature rules
include("D2Q4.jl")
include("D2Q5.jl")
include("D2Q9.jl")
include("D2Q13.jl")
include("D2Q17.jl")
include("D2Q21.jl")
include("D2Q37.jl")

const Quadratures = (
    D2Q4 = D2Q4(),
    D2Q5 = D2Q5(),
    D2Q9 = D2Q9(),
    D2Q13 = D2Q13(),
    D2Q17 = D2Q17(),
    D2Q21 = D2Q21(),
    D2Q37 = D2Q37(),
)

include("plot-quadratures.jl")
