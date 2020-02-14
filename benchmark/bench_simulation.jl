module BenchSimulation

using BenchmarkTools
using LatticeBoltzmann

import LatticeBoltzmann: collide!, stream!, apply_boundary_conditions!, next!, simulate

suite = BenchmarkGroup()

function initialize_model(q = D2Q9(), τ = 1.0, scale = 2)
    cs = 1 / q.speed_of_sound_squared
    benchmark_problem = LatticeBoltzmann.CouetteFlow(cs * (τ - 0.5), scale)
    return LatticeBoltzmann.LatticeBoltzmannModel(
        benchmark_problem,
        q,
        process_method = LatticeBoltzmann.ProcessingMethod(benchmark_problem, false, 0),
    )
end
for q in LatticeBoltzmann.Quadratures
    suite[string(q)] = BenchmarkGroup([string(q), "simulation"])

    model = initialize_model(q)

    suite[string(q)]["collision"] =
        @benchmarkable(collide!(model, time = 0.0), setup = (model = initialize_model($q)))
    suite[string(q)]["stream"] =
        @benchmarkable(stream!(model), setup = (model = initialize_model($q)))
    suite[string(q)]["boundary conditions"] = @benchmarkable(
        apply_boundary_conditions!(model, time = 0.0),
        setup = (model = initialize_model($q))
    )
    suite[string(q)]["processing"] =
        @benchmarkable(next!(model, 0), setup = (model = initialize_model($q)))
    suite[string(q)]["simulate"] =
        @benchmarkable(simulate(model, 0:10), setup = (model = initialize_model($q)))
end

end  # module
BenchSimulation.suite
