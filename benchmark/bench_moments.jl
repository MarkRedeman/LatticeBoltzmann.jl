module BenchMoments

using BenchmarkTools
using LatticeBoltzmann
using StaticArrays

import LatticeBoltzmann:
    density,
    velocity!,
    temperature,
    pressure,
    momentum_flux,
    deviatoric_tensor

suite = BenchmarkGroup()

function initialize_benchmark(q = D2Q9(), τ = 1.0, scale = 2)
    benchmark_problem = LatticeBoltzmann.TGV(q, τ, scale)

    LatticeBoltzmann.initialize(
        LatticeBoltzmann.ZeroVelocityInitialCondition(),
        q,
        benchmark_problem
    )
end

for q = LatticeBoltzmann.Quadratures
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

    suite[string(q)]["pressure"] =
        @benchmarkable(
            pressure($q, f, ρ, u),
            setup = (
                f = copy($q.weights);
                ρ = density($q, f);
                u = zeros(eltype(f), dimension($q));
                velocity!($q, f, ρ, u)
            )
        )

    suite[string(q)]["temperature"] =
        @benchmarkable(
            pressure($q, f, ρ, u),
            setup = (
                f = copy($q.weights);
                ρ = density($q, f);
                u = zeros(eltype(f), dimension($q));
                velocity!($q, f, ρ, u)
            )
        )

    suite[string(q)]["momentum_flux"] =
        @benchmarkable(
            momentum_flux($q, f, ρ, u),
            setup = (
                f = copy($q.weights);
                ρ = density($q, f);
                u = zeros(eltype(f), dimension($q));
                velocity!($q, f, ρ, u)
            )
        )

    suite[string(q)]["deviatoric_tensor"] =
        @benchmarkable(
            deviatoric_tensor($q, τ, f, ρ, u),
            setup = (
                τ = 1.0;
                f = copy($q.weights);
                ρ = density($q, f);
                u = zeros(eltype(f), dimension($q));
                velocity!($q, f, ρ, u)
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
