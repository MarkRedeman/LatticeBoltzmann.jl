using LatticeBoltzmann, Plots, LaTeXStrings, DataFrames, JLD2
import LatticeBoltzmann: StopCriteria,
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


line_style(q::D2Q4) = (:line)
line_style(q::D2Q5) = (:line, :dash, 1.0, 3.0)
line_style(q::D2Q9) = (:line)
line_style(q::D2Q13) = (:line)
line_style(q::D2Q17) = (:line)
line_style(q::D2Q21) = (:line)
line_style(q::D2Q37) = (:line, :dash, 1.0, 3.0)
marker_style(q::Quadrature) = ()

marker_style(q::D2Q4) = ()
marker_style(q::D2Q5) = (:hexagon)
marker_style(q::D2Q9) = ()
marker_style(q::D2Q13) = ()
marker_style(q::D2Q17) = ()
marker_style(q::D2Q21) = ()
marker_style(q::D2Q37) = (:hexagon)

function shear_wave_convergence_analysis(
    q = D2Q9(),
    initialization_strategy = AnalyticalEquilibrium(),
    τ_lb = 1.0;
    scales = [1, 2, 4, 8],
)
    ν_lb = τ_lb / (2.0 * q.speed_of_sound_squared)
    t_end = 1.0

    results = map(scales) do scale
        problem = DecayingShearFlow(ν_lb, scale, static = true)

        Δt = delta_t(problem)
        n_steps = round(Int, t_end / Δt)

        # CompareWithAnalyticalSolution
        process_method = TrackHydrodynamicErrors(
            problem,
            false,
            n_steps,
            LatticeBoltzmann.NoStoppingCriteria(),
        )

        @time res = simulate(
            problem,
            q,
            process_method = process_method,
            initialization_strategy = initialization_strategy,
            t_end = t_end,
        )

        return res.processing_method.df[end]
    end

    return (quadrature = q, scales = scales, results = results)
end
