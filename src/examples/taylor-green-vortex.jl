module Example
module TaylorGreenVortex

using BenchmarkTools
using StatsPlots
using DataFrames
using lbm
using Plots


results = let
    stats = DataFrame([Float64[], Int[], Any[]], [:nu, :scale, :stats])

    quadratures = [
        D2Q4(),
        # D2Q5(),
        D2Q9(),
        # D2Q17(),
    ]

    quadrature = last(quadratures)
    example = TaylorGreenVortexExample(1.0 / 6.0, 2)
    # example = DecayingShearFlow(1.0 / 6.0, 4)

    result = lbm.siumlate(example, quadrature)
    return result


    νs = (0.0:0.5:4.0) ./ 6.0
    scales = [1, 2, 4, 6, 8, 10, 12]

    νs = (0.0:0.5:4.0) ./ 6.0
    scales = [1, 2]
    for ν in νs
    for scale = scales
    # scale = 1
    # ν = 1.0 / 6.0
    # ν = 0.0
    # scale = 2
    # ν = 1.0 / 6.0
    example = TaylorGreenVortexExample(ν, scale)

    result = lbm.siumlate(example, quadrature);
    push!(stats, [ν, scale, result[2]])
    end
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
    plot!(
        scale=:log10,
        legend=:bottomleft,
        legendfontsize=5
    )
    gui()

    # nu_idx = 0
    # plot(nu_scale_error[nu_idx * length(scales) .+ (1:length(scales)), 2], nu_scale_error[nu_idx * length(scales) .+ (1:length(scales)), 3])
    return stats
end

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
