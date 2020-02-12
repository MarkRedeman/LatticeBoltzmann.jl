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


function initial_conditions_analysis()
    q = D2Q9()
    scales = [1, 2, 4]

    Nx = 31
    Ny = 17
    τs = [0.51, 0.6, 0.8, 1.0]

    scale = 2

    Nx = 31
    Ny = 17

    τ = 0.8
    t_d = 840

    iteration_strategies = [
        LatticeBoltzmann.AnalyticalVelocityAndStress(),
        LatticeBoltzmann.AnalyticalEquilibriumAndOffEquilibrium(),
    ]
    iteration_strategies = [
        # With this initialization strategy it is assumed that the initial density field
        # (or equivalently the pressure field p) is not available
        LatticeBoltzmann.ConstantDensity(),

        # TODO
        LatticeBoltzmann.AnalyticalVelocityAndStress(),

        # The initial pressure field p is known, which, using the equation fo state,
        # is used to set the initial density
        LatticeBoltzmann.AnalyticalEquilibrium(),

        # Both an initial pressure field and gradient of the velocity is known.
        # The gradient of the velocity is used to initialize the offequilibrium components
        LatticeBoltzmann.AnalyticalEquilibriumAndOffEquilibrium(),

        # We initialize f using an iterative procedure where only the density is conserved
        # it was shown that this procedure gives consistent itnitial conditions for both
        # the equilibrium and off equilibrium components
        LatticeBoltzmann.IterativeInitializationMeiEtAl(τ, 1E-10),
        LatticeBoltzmann.IterativeInitializationMeiEtAl(1.0, 1E-10),
    ]

    init_res = map(iteration_strategies) do initialization
        u_0 = sqrt(0.01)
        u_0 = 0.03
        Nx = 96
        Ny = 72
        τ = 0.8

        # u_0 = 0.06
        # Nx = 48
        # Ny = 36

        # u_0 = 0.01
        # Nx = 32
        # Ny = 32
        # ν = 0.002
        # τ = q.speed_of_sound_squared * ν + 0.5

        problem = LatticeBoltzmann.TGV(q, τ, scale, Nx, Ny, u_0)
        t_end = round(Int, LatticeBoltzmann.decay_time(problem))
        @show t_end

        @time model = LatticeBoltzmann.LatticeBoltzmannModel(
            problem,
            q,
            initialization_strategy = initialization,
            process_method = LatticeBoltzmann.ProcessingMethod(problem, true, t_end),
        )

        @time solution = simulate(model, 1:t_end)
    end

    return init_res
end
