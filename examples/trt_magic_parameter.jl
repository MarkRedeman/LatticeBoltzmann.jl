module TRT_MAGIC_EXAMPLE

using LatticeBoltzmann, Plots, LaTeXStrings, DataFrames, JLD2
using Logging, TerminalLoggers, ProgressMeter
import LatticeBoltzmann:
    StopCriteria,
    CompareWithAnalyticalSolution,
    TrackHydrodynamicErrors,
    ZeroVelocityInitialCondition,
    IterativeInitializationMeiEtAl,
    ConstantDensity,
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
        LatticeBoltzmann.VelocityConvergenceStoppingCriteria(1E-7, problem),
    )

    collision_model = LatticeBoltzmann.TRT(
        τ_s,
        τ_a,
        (x_idx, y_idx, t) -> LatticeBoltzmann.lattice_force(problem, x_idx, y_idx, t),
    )

    result = LatticeBoltzmann.simulate(
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

function main(
    quadratures = [D2Q9()],
    scale = 2,
    τ_s_range = range(0.51, stop = 1.0, step = 0.1),
    τ_a_range = τ_s_range,
)
    with_logger(TerminalLogger(stderr, Logging.Warn)) do
        s = []
        for q in quadratures
            @showprogress "Computing optimal relaxation times..." for τ_s in τ_s_range,
                τ_a in τ_a_range

                ν = (τ_s - 0.5) / q.speed_of_sound_squared
                problem = PoiseuilleFlow(ν, scale)

                result = solve(problem, q, τ_s, τ_a)
                errors = result.processing_method.df[end]

                push!(
                    s,
                    (
                        τ_s = τ_s,
                        τ_a = τ_a,
                        quadrature = q,
                        error_u = errors.error_u,
                        error_p = errors.error_p,
                        error_σ_xx = errors.error_σ_xx,
                        error_σ_xy = errors.error_σ_xy,
                    ),
                )
            end
        end
        return s
    end
end
function plot_results(x)
    τ_as = map(d -> d.τ_a, x) |> unique
    τ_ss = map(d -> d.τ_s, x) |> unique
    p = contour(
        τ_ss,
        τ_as,
        (s, a) -> x[findfirst(d -> d.τ_s == s && d.τ_a == a, x)].error_u,
        xlabel = L"\tau_s",
        ylabel = L"\tau_a",
        fill = true,
    )
    return p

    # @pgf Axis(
    #     {
    #         view = (0, 90),
    #         colorbar,
    #         "colormap/jet"
    #     },
    #     Plot3(
    #         {
    #             surf,
    #             shader = "flat",

    #         },
    #         Table(getfield.(x, :τ_s), getfield.(x, :τ_a), getfield.(x, :error_u)))
    # )
end

end
