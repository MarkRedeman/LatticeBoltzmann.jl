module BenchMoments

using BenchmarkTools
using LatticeBoltzmann
using StaticArrays

import LatticeBoltzmann:
    density,
    velocity!

suite = BenchmarkGroup()

function initialize_benchmark(q = D2Q9(), τ = 1.0, scale = 2)
    benchmark_problem = LatticeBoltzmann.TGV(q, τ, scale)

    LatticeBoltzmann.initialize(
        LatticeBoltzmann.ZeroVelocityInitialCondition(),
        q,
        benchmark_problem
    )
end

for q = [D2Q9(), D2Q13(), D2Q17(), D2Q21(), D2Q37()]
    suite[string(q)] = BenchmarkGroup([string(q), "moments"])

    suite[string(q)]["density"] =
        @benchmarkable(
            density($q, f),
            setup = (f = copy($q.weights);)
        )

    suite[string(q)]["velocity"] =
        @benchmarkable(
            velocity!($q, f, ρ, u),
            setup = (
                f = copy($q.weights);
                ρ = density($q, f);
                u = zeros(eltype(f), dimension($q))
            )
        )


    suite[string(q)]["density", "static"] =
        @benchmarkable(
            density(static_q, f),
            setup = (
                f = $q.weights |> MVector{length($q.weights)};
                static_q = D2Q9(
                    $q.abscissae |> SMatrix{dimension($q), length($q.weights)},
                    tuple($q.weights...) |> SVector{length($q.weights)},
                    $q.speed_of_sound_squared
                )
            )
        )
    suite[string(q)]["velocity", "static"] =
        @benchmarkable(
            velocity!(static_q, f, ρ, u),
            setup = (
                f = $q.weights |> MVector{length($q.weights)};
                ρ = density($q, f);
                u = @MVector zeros(eltype(f), dimension($q));
                static_q = D2Q9(
                    $q.abscissae |> SMatrix{dimension($q), length($q.weights)},
                    tuple($q.weights...) |> SVector{length($q.weights)},
                    $q.speed_of_sound_squared
                )
            )
        )
end

end  # module
BenchMoments.suite
