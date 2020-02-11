using LatticeBoltzmann, Plots, DataFrames

# LatticeBoltzmann.analyze_convergence(D2Q9(), (scale, viscosity) -> TaylorGreenVortex(viscosity, scale, static = true), 1.0 / 6.0, 2)
# LatticeBoltzmann.analyze_convergence(D2Q9(), (scale, viscosity) -> TaylorGreenVortex(viscosity, scale, static = false), 1.0 / 6.0, 3)
# LatticeBoltzmann.analyze_convergence(D2Q9(), (scale, viscosity) -> PoiseuilleFlow(viscosity, scale, static = true), 1.0 / 6.0, 3)
# LatticeBoltzmann.analyze_convergence(D2Q9(), (scale, viscosity) -> DecayingShearFlow(viscosity, scale, static = true), 1.0 / 6.0, 3)

stats = DataFrame([Float64[], Int[], Any[]], [:nu, :scale, :stats])

quadratures = LatticeBoltzmann.Quadratures

# ProcessingMethod(problem::DecayingShearFlow, should_process, n_steps, stop_criteria = StopCriteria(problem)) =
#     CompareWithAnalyticalSolution(problem, should_process, n_steps, stop_criteria)
let
    q = D2Q9()
    τ = 1.0 / 6.0
    scale = 1

    problem = DecayingShearFlow(τ, scale, static = true)
    result = LatticeBoltzmann.simulate(problem, q)

    problem = DecayingShearFlow(τ, scale, static = false)
    result = LatticeBoltzmann.simulate(problem, q)
end

scales = 1:8
x = map(
    scale -> begin
        result = LatticeBoltzmann.simulate(
            DecayingShearFlow(τ, scale, static = true),
            D2Q9(),
            t_end = 0.5,
            should_process = false,
        )
        return result.processing_method.df[end]
    end,
    scales,
)

plot(
    plot(
        map(scale -> 16 * scale, scales),
        getfield.(x, :error_u),
        scale = :log10,
        title = "u",
    ),
    plot(
        map(scale -> 16 * scale, scales),
        getfield.(x, :error_σ_xx),
        scale = :log10,
        title = "sigma xx",
    ),
    plot(
        map(scale -> 16 * scale, scales),
        getfield.(x, :error_σ_xy),
        scale = :log10,
        title = "sigma xy",
    ),
    plot(
        map(scale -> 16 * scale, scales),
        getfield.(x, :error_p),
        scale = :log10,
        title = "p",
    ),
)
