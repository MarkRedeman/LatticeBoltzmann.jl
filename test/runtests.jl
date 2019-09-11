using Test

@test 1 == 1

using lbm

quadratures = [
    D2Q4(),
    D2Q5(),
    D2Q9()
]

@testset "moments of the equilibrium distribution of the $q lattice" for q in quadratures
    ρ = 1
    u = [1 1]
    T = 1
    @test 1 == 1
end
