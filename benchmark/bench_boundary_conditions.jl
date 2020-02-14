module BenchBoundaryConditions

using BenchmarkTools
using LatticeBoltzmann

import LatticeBoltzmann: apply!, Direction, North, East, South, West, BounceBack, MovingWall

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
    suite[string(q)] = BenchmarkGroup([string(q), "boundary conditions"])

    τ = 1.0
    f = initialize_benchmark(q, τ)
    NX, NY, NF = size(f)
    for direction in [North(), East(), South(), West()]
        suite[string(q)]["bounce back", direction] = @benchmarkable(
            apply!(bc, $q, f_new, f_old),
            setup = (bc = BounceBack($direction, 1:($NX), 1:($NY));
            f_new = copy($f);
            f_old = copy($f))
        )
    end
end

end  # module
BenchBoundaryConditions.suite
