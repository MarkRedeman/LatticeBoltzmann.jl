using lbm, Plots, DataFrames

# lbm.analyze_convergence(D2Q9(), (scale, viscosity) -> TaylorGreenVortexExample(viscosity, scale, static = true), 1.0 / 6.0, 2)
# lbm.analyze_convergence(D2Q9(), (scale, viscosity) -> TaylorGreenVortexExample(viscosity, scale, static = false), 1.0 / 6.0, 3)
# lbm.analyze_convergence(D2Q9(), (scale, viscosity) -> PoiseuilleFlow(viscosity, scale, static = true), 1.0 / 6.0, 3)
# lbm.analyze_convergence(D2Q9(), (scale, viscosity) -> DecayingShearFlow(viscosity, scale, static = true), 1.0 / 6.0, 3)

stats = DataFrame([Float64[], Int[], Any[]], [:nu, :scale, :stats])

quadratures = lbm.Quadratures

let
    q = D2Q9()
    τ = 1.0 / 6.0
    scale = 1

    problem = DecayingShearFlow(τ, scale, static = true)
    result = lbm.siumlate(problem, q)

    problem = DecayingShearFlow(τ, scale, static = false)
    result = lbm.siumlate(problem, q)
end

scales = 1:8
x = map(
    scale -> begin
    result = lbm.siumlate(
        DecayingShearFlow(τ, scale, static = true),
        D2Q9(),
        t_end = .5,
        should_process=false
    );
    return result.processing_method.df[end]
    end,
    scales
)

plot(
plot(map(scale -> 16*scale, scales), getfield.(x, :error_u), scale=:log10, title="u"),
plot(map(scale -> 16*scale, scales), getfield.(x, :error_σ_xx), scale=:log10, title="sigma xx"),
plot(map(scale -> 16*scale, scales), getfield.(x, :error_σ_xy), scale=:log10, title="sigma xy"),
plot(map(scale -> 16*scale, scales), getfield.(x, :error_p), scale=:log10, title="p"),
)
