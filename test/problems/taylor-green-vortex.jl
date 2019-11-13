@testset "Taylor Green Vortex Problem" begin #for q in quadratures
    q = D2Q9()

    scale = 1
    ν = 1.0 / 6.0
    problem = TaylorGreenVortexExample(ν, scale,)
    f = lbm.initialize(q, problem)

    ρ = lbm.density(q, f)
    j = lbm.momentum(q, f)
    T = lbm.temperature(q, f, ρ, j ./ ρ)

    # TODO: when an a thermal model is considered, the density fluctuates
    # TODO: when a thermal model is considerd we have constant density
end

# struct LatticeNode{Q, T}
#     distributions::Array{Q.d}
#     ρ::Float64
#     u::Vector{Float64}
#     T::Float64


# end
