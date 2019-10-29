module Example
module PoiseuilleFlowExample

using BenchmarkTools
using StatsPlots
using DataFrames
using lbm
using Plots


quadratures = [
    D2Q4(),
    # D2Q5(),
    D2Q9(),
    D2Q17(),
]

q = last(quadratures)
scale = 1

# relaxation time (BGK model) (tau=sqrt(3/16)+0.5 gives exact solution)
# tau=sqrt(3/16)+0.5;
# nu=(2*tau-1)/6;

let
    # ν = (1.4383 - 0.5) / q.speed_of_sound_squared
    q = D2Q9()
    stats = DataFrame([Float64[], Float64[]], [:τ, :u_error])
    for τ = 0.51:0.01:2.0
        ν = (τ - 0.5) / q.speed_of_sound_squared

        problem = PoiseuilleFlow(ν, 1, static = true)

        result = lbm.siumlate(
            problem,
            q,
            t_end = 100.0,
            should_process = false
        )
        if (! isnan(result[2].u_error[end]))
            push!(stats, [τ, result[2].u_error[end]])
        end
    end

    plot(stats.τ, stats.u_error, yscale=:log10)
    gui()

    @show stats.τ[argmin(stats.u_error)]
# stats.τ[argmin(stats.u_error)] = 1.438
# 1.438
# julia> @show stats.τ[argmin(stats.u_error)]
# stats.τ[argmin(stats.u_error)] = 1.76
# 1.76
end
# @show result[2]
    # return result


    νs = (0.0:2.0:6.0) ./ 6.0
    # νs = [1.0 / 6.0]
    scales = [1, 2, 4, 6, 8, 10, 12]#, 8, 10, 12]
    scales = [1, 2, 4] #, 8, 16]#, 8, 10, 12]

    # νs = (0.0:0.5:4.0) ./ 6.0
    # scales = [1, 2, 4]
    # νs = [1.0 / 6.0]
    for ν in νs, scale in scales
        continue

        # example = TaylorGreenVortexExample(ν, scale, static = false)
        example = DecayingShearFlow(ν, scale, static = true)
        # example = TaylorGreenVortexExample(τ, scale, static = true)

        result = lbm.siumlate(example, quadrature, base = 20);
        # @show result[2]
        push!(stats, [ν, scale, result[2]])
    end
    # return stats

    # end
    # end

    s = stats
    using Plots, LaTeXStrings
    nu_scale_error = Array{Float64}([s.nu s.scale map(i -> i.u_error[end], s.stats)])
    plot()
    for nu_idx = 1:length(νs)
        nu = nu_scale_error[nu_idx * length(scales), 1]
        nu_round = round(nu, digits = 2)

        plot!(
            map(
                i -> 16 * i,
                nu_scale_error[(nu_idx - 1) * length(scales) .+ (1:length(scales)), 2]
            ),
            nu_scale_error[(nu_idx - 1) * length(scales) .+ (1:length(scales)), 3],
            label="nu = $nu_round"
        )
    end
    plot!(x -> 10.0 * x^(-2), label="y(x) = 10x^(-2)", linestyle=:dash)
    # plot!(x -> 10.0 * x^(-3), label="y(x) = 10x^(-3)", linestyle=:dash)
    # plot!(x -> 10.0 * x^(-4), label="y(x) = 10x^(-4)", linestyle=:dash)
    # plot!(x -> 10.0 * x^(-5), label="y(x) = 10x^(-5)", linestyle=:dash)
    plot!(
        scale=:log10,
        legend=:bottomleft,
        legendfontsize=5
    )
    gui()

    @show -log.(nu_scale_error[2:end, 3] ./ nu_scale_error[1:(end - 1), 3]) ./ log.(nu_scale_error[2:end, 2] ./ nu_scale_error[1:(end - 1), 2])


    # nu_idx = 0
    # plot(nu_scale_error[nu_idx * length(scales) .+ (1:length(scales)), 2], nu_scale_error[nu_idx * length(scales) .+ (1:length(scales)), 3])

# Some thoughs: hide the storage of the distributions f inside of an interface
# so that we can do:
# for x in xs, y in ys, z in zs
#  f = fs[x, y ,z]::Vector
#
# the interface could hide it by returning
# return @view f_internal[x, y, z, :]


# Idea: introduce an Initial Value Problem
# problem = InitialValueProblem
# solution = solve(problem, LBM(Lattice, CollisionModel))
# LBM(Lattice, CollisionModel) can be a solution method

end
end
