import Base: OneTo
import LatticeBoltzmann: apply!, Direction, North, East, South, West, BounceBack, MovingWall

function initialize_benchmark(q = D2Q9(), τ = 1.0, scale = 2)
    benchmark_problem = LatticeBoltzmann.TGV(q, τ, scale)

    LatticeBoltzmann.initialize(LatticeBoltzmann.ZeroVelocityInitialCondition(), q, benchmark_problem)
end

@testset "Boundary conditions with $q" for q in [D2Q9(), D2Q13(), D2Q17(), D2Q21(), D2Q37()]
    τ = 1.0
    f = initialize_benchmark(q, τ)
    NX, NY, NF = size(f)
    for direction = [North(), East(), South(), West()]
        bc = BounceBack(direction, 1:NX, 1:NY)
        f_new = copy(f)
        f_old = copy(f)

        apply!(bc, q, f_new, f_old)

        # Assert that the tangential velocity is 0
    end
end
