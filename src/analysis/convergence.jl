using StatsPlots
using DataFrames
using lbm
using Plots

function test_convergence(problem::Problem, scales) where { Problem <: InitialValueProblem }
end

function analyze_convergence(
    q::lbm.Quadrature,
    p,
    viscosity::Float64,
    N = 2,
    t_end = 2pi
)
    stats = DataFrame([Float64[], Int[], Float64[]], [:nu, :scale, :u_error])

    for scale = 0:N
        problem = p(2^scale, viscosity)
        result = lbm.siumlate(problem, q, base = 200, should_process = false, t_end = t_end)
        push!(stats, [viscosity, 16 .* 2^scale, result[2].u_error[end]])
    end

    plot_convergence(stats, viscosity)
    gui()

    @show -log.(
        stats.u_error[2:end] ./ stats.u_error[1:(end - 1)]
    ) ./ log.(
        stats.nu[2:end] ./ stats.nu[1:(end - 1)]
    )

    return stats
end
function plot_convergence(stats, viscosity)
    nu_round = round(viscosity, digits = 2)

    p = plot()
    plot!(p, stats.scale, stats.u_error , label="nu = $nu_round")
    plot!(p, x -> stats.u_error[1] * 10.0 * x^(-2), label="y(x) = 10x^(-2)", linestyle=:dash)
    plot!(p, scale=:log10, legend=:bottomleft, legendfontsize=5)
    return p
end
