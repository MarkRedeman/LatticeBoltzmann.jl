using lbm, Plots, DataFrames

function plot_convergence_analysis_of(x)
p = plot()
plot!(p, map(scale -> 16*scale, scales), getfield.(x, :error_u), scale=:log10, label="u")
plot!(p, map(scale -> 16*scale, scales), getfield.(x, :error_σ_xx), scale=:log10, label="sigma xx")
plot!(p, map(scale -> 16*scale, scales), getfield.(x, :error_σ_xy), scale=:log10, label="sigma xy")
plot!(p, map(scale -> 16*scale, scales), getfield.(x, :error_p), scale=:log10, label="p")
plot!(p, map(scale -> 16*scale, scales), x -> 0.1 * x.^(-2), scale=:log10, label="x^-2")
plot!(p, legend=:bottomleft)
    p
end
function plot_convergence_rate!(rate = -2, p = nothing)
    plot!(x -> 0.1 * x.^(rate), label="x^" * string(rate))
end

q = D2Q9(); @time x = map(
    scale -> begin
    result = lbm.siumlate(
        DecayingShearFlow(1.000/ (q.speed_of_sound_squared * 2), scale, static = false),
        q,
        t_end = 0.1250,
        should_process=false
    );
    return result.processing_method.df[end]
    end,
    scales
) |> plot_convergence_analysis_of

q = D2Q9();
scale = 2//1;
τ = 1.00 / (2 *q.speed_of_sound_squared);
result = lbm.siumlate(TaylorGreenVortex(τ, scale, static = false), q, t_end = 0.1250, should_process=true);

Nx = 31
Ny = 17
τs = [0.51, 0.6, 0.8, 1.0]
τ = last(τs)

q = D2Q9()
scales = [1, 2, 4, 8, 16]
ρ_0 = 1.0

A = 1.0
B = -1.0
a = 1.0
b = 1.0

for scale in scales
    u_max = √(0.001) / scale
    ν = τ / (2.0 * q.speed_of_sound_squared)
    NX = Int(scale * 31)
    NY = Int(scale * 17)

    problem = TaylorGreenVortex(ρ_0, u_max, ν, NX, NY, (2pi, 2pi), false, A, B, a, b)

    result = lbm.siumlate(
        problem,
        q,
        t_end = 1.0,
        should_process=false
    )
end

scales = [1//2, 1, 2, 4, 8]
@time x = map(
    scale -> begin
    u_max = √(0.0001) / scale
    ν = τ / (2.0 * q.speed_of_sound_squared)
    NX = Int(scale * 16)
    NY = Int(scale * 16)

    problem = TaylorGreenVortex(ρ_0, u_max, ν, NX, NY, (2pi, 2pi), false, A, B, a, b)

    result = lbm.siumlate(
        problem,
        q,
        # t_end = 1.0,
        t_end = 0.1250,
        should_process=false
    )
    return result.processing_method.df[end]
    end,
    scales
) |> plot_convergence_analysis_of

q = D2Q9(); @time x = map(
           scale -> begin
           result = lbm.siumlate(
               TaylorGreenVortex(1.000/ (q.speed_of_sound_squared * 2), scale, static = false),
               q,
               t_end = 0.1250,
               should_process=false
           );
           return result.processing_method.df[end]
           end,
           scales
       ) |> plot_convergence_analysis_of
