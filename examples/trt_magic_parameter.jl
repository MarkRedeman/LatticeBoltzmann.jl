module TRT_MAGIC_EXAMPLE

using LaTeXStrings
using DataFrames
using LatticeBoltzmann
using Plots
using JLD2

import LatticeBoltzmann: StopCriteria,
    CompareWithAnalyticalSolution,
    TrackHydrodynamicErrors,
    ZeroVelocityInitialCondition,
    IterativeInitializationMeiEtAl,
    density,
    velocity!,
    dimensionless_velocity,
    ProcessingMethod,
    next!,
    InitializationStrategy,
    ShowVelocityError

# Choose an end time such that we can be certain that a steady state solution
# has been found before terminating
const T_END = 100.0

function solve(problem, q, τ_s, τ_a)
    t_end = T_END
    Δt = LatticeBoltzmann.delta_t(problem)
    n_steps = round(Int, t_end / Δt)

    process_method = TrackHydrodynamicErrors(
        problem,
        false,
        n_steps,
        LatticeBoltzmann.VelocityConvergenceStoppingCriteria(1E-7, problem)
    )

    collision_model = LatticeBoltzmann.TRT(
        τ_s,
        τ_a,f
        (x_idx, y_idx, t) -> LatticeBoltzmann.lattice_force(problem, x_idx, y_idx, t)
    )

    @time result = LatticeBoltzmann.simulate(
        problem,
        q,
        t_end = t_end,
        should_process = false,
        collision_model = collision_model,
        process_method = process_method,
        initialization_strategy = ZeroVelocityInitialCondition(),
    )

    return result
end

function plot_main(result; Quadratures = [D2Q9()], scale = 2)
    s = []
    for q in Quadratures
        for τ_s in range(0.51, stop = 10.0, step = 0.05),
            τ_a in range(0.51, stop = 10.0, step = 0.05)

            ν = (τ_s - 0.5) / q.speed_of_sound_squared
            problem = PoiseuilleFlow(ν, scale)

            result = solve(problem, q, τ_s, τ_a)
            errors = result.processing_method.df[end]

            push!(s, (
                τ_s = τ_s,
                τ_a = τ_a,
                quadrature = q,
                error_u = errors.error_u,
                error_p = errors.error_p,
                error_σ_xx = errors.error_σ_xx,
                error_σ_xy = errors.error_σ_xy
            ))
        end
    end
    return s
end

function plot_results(x)
    τ_as = map(d -> d.τ_a, x) |> unique
    τ_ss = map(d -> d.τ_s, x) |> unique
    p = contour(
        τ_ss,
        τ_as,
        (s, a) -> x[findfirst(d -> d.τ_s == s && d.τ_a == a, x)].error_u,
        xlabel=L"\tau_s",
        ylabel=L"\tau_a",
        fill=true
    )
    return p
end

end
