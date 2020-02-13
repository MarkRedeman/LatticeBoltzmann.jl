module BenchCollisionModels

using BenchmarkTools
using LatticeBoltzmann

suite = BenchmarkGroup()

function initialize_benchmark(q = D2Q9(), τ = 1.0, scale = 2)
    benchmark_problem = LatticeBoltzmann.TGV(q, τ, scale)

    LatticeBoltzmann.initialize(LatticeBoltzmann.ZeroVelocityInitialCondition(), q, benchmark_problem)
end

for q = [D2Q9(), D2Q13(), D2Q17(), D2Q21(), D2Q37()]
    suite[string(q)] = BenchmarkGroup([string(q), "collision_models"])

    srt_equilibrium = SRT(1.0)
    trt_equilibrium = TRT(1.0, 1.0)
    mrt_equilibrium = MRT(q, 1.0)

    τ = 1.0
    f = initialize_benchmark(q, τ)

    suite[string(q)]["srt", 1.0] =
        @benchmarkable(
            collide!(cm, $q, f_in, f_out),
            setup = (cm = $srt_equilibrium; f_in = copy($f); f_out = copy($f))
        )

    suite[string(q)]["srt", 1.3] =
        @benchmarkable(
            collide!(cm, $q, f_in, f_out),
            setup = (cm = SRT(1.3); f_in = copy($f); f_out = copy($f))
        )

    suite[string(q)]["srt", 0.8] =
        @benchmarkable(
            collide!(cm, $q, f_in, f_out),
            setup = (cm = SRT(0.8); f_in = copy($f); f_out = copy($f))
        )

    suite[string(q)]["trt", (1.0, 1.0)] =
        @benchmarkable(
            collide!(cm, $q, f_in, f_out),
            setup = (cm = $trt_equilibrium; f_in = copy($f); f_out = copy($f))
        )

    suite[string(q)]["trt", (0.9, 1.1)] =
        @benchmarkable(
            collide!(cm, $q, f_in, f_out),
            setup = (cm = TRT(0.9, 1.1); f_in = copy($f); f_out = copy($f))
        )

    suite[string(q)]["mrt-equilibrium", 1.0] =
        @benchmarkable(
            collide!(cm, $q, f_in, f_out),
            setup = (cm = $mrt_equilibrium; f_in = copy($f); f_out = copy($f))
        )

    scale = 2
    benchmark_problem = LatticeBoltzmann.DecayingShearFlow(1.0 / (2 * q.speed_of_sound_squared), scale, static = true)
    f_force = LatticeBoltzmann.initialize(LatticeBoltzmann.ZeroVelocityInitialCondition(), q, benchmark_problem)
    suite[string(q)]["srt", "force"] =
        @benchmarkable(
            collide!(cm, $q, f_in, f_out),
            setup = (
                cm = CollisionModel(SRT, $q, $benchmark_problem);
                f_in = copy($f_force);
                f_out = copy($f_force)
            )
        )

    suite[string(q)]["trt", "force"] =
        @benchmarkable(
            collide!(cm, $q, f_in, f_out),
            setup = (
                cm = CollisionModel(TRT, $q, $benchmark_problem);
                f_in = copy($f_force);
                f_out = copy($f_force)
            )
        )
end

end  # module
BenchCollisionModels.suite
