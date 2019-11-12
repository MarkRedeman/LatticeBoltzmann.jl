using LinearAlgebra

abstract type Quadrature end
# abstract type Quadrature{Thermal} end
abstract type Lattice end

DdQq(d, q) = "HOI"

"""
Get the index of the abscissae pointing in the opposite direction of the given index

Here we assume that the abscissae of a quadrature are ordered such that each abscissae
follows its opposite abscissae
"""
function opposite(q::Quadrature, idx::Int64)
    if idx == 1
        return 1
    end
    if (mod(idx, 2) == 0)
        return idx + 1
    end
    return idx - 1
end


δ(α, β) = α == β ? 1 : 0
include("hermite.jl")
include("equilibrium-coefficients.jl")
include("thermodynamics.jl")
include("equilibrium.jl")

# Specific quadrature rules
include("D2Q4.jl")
include("D2Q5.jl")
include("D2Q9.jl")
include("D2Q13.jl")
include("D2Q17.jl")
include("D2Q21.jl")
include("D2Q37.jl")


# plot(map(q -> lbm.plot_quadrature(q, false), [D2Q4(), D2Q5(), D2Q9(), D2Q13(), D2Q17(), D2Q21(), D2Q37()])...)

"""
Show the arrows of a quadrature, useful for debugging
"""
function plot_quadrature(q::Quadrature, color_each_group_separatly = false)
    unique_weights = unique(q.weights)
    colors = distinguishable_colors(length(unique_weights))

    limit = max(q.abscissae[1, :]...) + 1

    p = plot(
        title=string(q),
        legend=nothing,
        framestyle=:grid,
        tickfontcolor=:gray,
        xlims=(-limit, limit),
        ylims=(-limit, limit),
        gridstyle=:dash,
    )

    for f_idx = 1 : length(q.weights)
        color_idx = findfirst(isequal(q.weights[f_idx]), unique_weights)

        if q.abscissae[:, f_idx] != [0, 0]
            quiver!(
                p,
                [0.0], [0.0],
                quiver=(
                    [q.abscissae[1, f_idx]],
                    [q.abscissae[2, f_idx]]
                ),
                color=color_each_group_separatly ? colors[color_idx] : :black
            )
        end
    end

    for f_idx = 1 : length(q.weights)
        if q.abscissae[:, f_idx] != [0, 0]
            continue
        end
        color_idx = findfirst(isequal(q.weights[f_idx]), unique_weights)
        scatter!(
            p,
            [(q.abscissae[1, f_idx], q.abscissae[2, f_idx])],
            # color=colors[color_idx],
            color=color_each_group_separatly ? colors[color_idx] : :gray,
            foreground_color_border=:black,
            markersize=2
        )
    end

    return p
end
