module BenchEquilibria

using BenchmarkTools
using LatticeBoltzmann

import LatticeBoltzmann: equilibrium, equilibrium!

suite = BenchmarkGroup()

function initialize_benchmark(q = D2Q9(), τ = 1.0, scale = 2)
    benchmark_problem = LatticeBoltzmann.TGV(q, τ, scale)

    LatticeBoltzmann.initialize(
        LatticeBoltzmann.ZeroVelocityInitialCondition(),
        q,
        benchmark_problem,
    )
end

for q in LatticeBoltzmann.Quadratures
    suite[string(q)] = BenchmarkGroup([string(q), "maxwell boltzmann equilibrium"])

    suite[string(q)]["initial equilibrium"] = @benchmarkable(
        equilibrium($q, ρ, u, T),
        setup = (ρ = 1.0;
        u = zeros(eltype($q.weights), dimension($q));
        T = 1.0)
    )

    suite[string(q)]["compute equilibrium"] = @benchmarkable(
        equilibrium!($q, ρ, u, T, f),
        setup = (ρ = 1.0;
        u = zeros(eltype($q.weights), dimension($q));
        T = 1.0;
        f = zeros(eltype($q.weights), length($q.weights)))
    )
end

end  # module
BenchEquilibria.suite
