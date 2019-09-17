
density(q::Quadrature, f::Array{Float64, N}) where N = sum(f, dims=N)
density(q::Quadrature, f::Array{Float64, 1}) = sum(f)

function momentum(q::Quadrature, f::Array{Float64, 3})
    j = Array{Float64}(undef, size(f, 1), size(f, 2), dimension(q))

    @inbounds for x = 1 : size(f, 1), y = 1 : size(f, 2)
        j[x, y, :] = momentum(q, f[x, y, :])
    end

    return j
end
function momentum(q::Quadrature, f::Array{Float64, 1})::Array{Float64, 1}
    j = zeros(dimension(q))
    momentum!(q, f, j)
    return j
end

function momentum!(q::Quadrature, f::Array{Float64, 1}, j::Array{Float64, 1})
    @inbounds for d in 1 : dimension(q)
        j[d] = 0.0
        for idx = 1:length(f)
            j[d] += f[idx] * q.abscissae[d, idx]
        end
    end
    return
end

function momentum2!(q::Quadrature, f::Array{Float64, 1}, j::Array{Float64, 1})
    @inbounds for d in 1 : dimension(q)
        j[d] = dot(f, view(q.abscissae, d, :))
        # j[d] = f' * view(q.abscissae, d, :)
    end
    return
end

function total_energy(q::Quadrature, f)
    return sum([
        f[x, y, f_idx] * (q.abscissae[:, f_idx]' * q.abscissae[:, f_idx]) for x = 1 : size(f, 1), y = 1 : size(f, 2), f_idx = 1 : size(f, 3)
    ], dims = 3)
end

function internal_energy(q::Quadrature, f, ρ, u)
    return total_energy(q, f) - kinetic_energy(q, f, ρ, u)
end

function kinetic_energy(q::Quadrature, f, ρ, u)
    return [
        ρ[x, y] * (u[x, y, 1].^2 + u[x, y, 2].^2) for x = 1 : size(f, 1), y = 1 : size(f, 2)
    ]
end

function dimension(q::Quadrature)::Int
    return size(q.abscissae, 1)
end

function temperature(q::Quadrature, f, ρ, u)
    return internal_energy(q, f, ρ, u) * (2 / dimension(q))
end
