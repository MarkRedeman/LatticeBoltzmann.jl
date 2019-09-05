
abstract type Quadrature end
abstract type Lattice end

include("D2Q4.jl")
include("D2Q5.jl")
include("D2Q9.jl")
include("D2Q17.jl")

DdQq(d, q) = "HOI"

abstract type InitialValueProblem end
function initialize(quadrature::Quadrature, lattice::Lattice, problem::InitialValueProblem)

end

density(q::Quadrature, f) = sum(f, dims=3)
function momentum(q::Quadrature, f)
    j = zeros(size(f, 1), size(f, 2), dimension(q))
    for x = 1 : size(f, 1), y = 1 : size(f, 2), d = 1 : dimension(q)
        j[x, y, d] = sum([
            f[x, y, f_idx] * q.abscissae[d, f_idx] for f_idx = 1 : size(f, 3)
        ])
    end
    return j
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

function dimension(q::Quadrature)
    return size(q.abscissae, 1)
end

function temperature(q::Quadrature, f, ρ, u)
    return internal_energy(q, f, ρ, u) * (2 / dimension(q))
end

function equilibrium(q::Quadrature, f_in)
    # Density
    f_ρ = density(q, f_in)

    # Momentum
    j = momentum(q, f_in)

    T = 1.0 #* temperature(quadrature, f_in, f_ρ, j ./ f_ρ)
    return equilibrium(q, f_ρ, j ./ f_ρ, T);
end

function equilibrium(q::Quadrature, ρ, u, T)
    f = zeros(size(ρ,1), size(ρ,2), length(q.weights));

    u_squared = u[:, :, 1].^2 + u[:, :, 2].^2
    for idx = 1:length(q.weights)
        cs = q.abscissae[1, idx] .* u[:, :, 1] .+ q.abscissae[2, idx] .* u[:, :, 2]

        f[:, :, idx] = _equilibrium(q, ρ, q.weights[idx], cs, u_squared, T, q.abscissae[1, idx]^2 + q.abscissae[2, idx]^2)
    end

    return f
end


function _equilibrium(q, ρ, weight, u_dot_xi, u_squared, T, xi_squared)
    # Truncated upto order 2
    sss = q.speed_of_sound_squared
    # return ρ .* weight .* (1.0 .+ sss * u_dot_xi .+ 4.5 * (u_dot_xi .* u_dot_xi) .- 1.5 * u_squared)
    # a = (T .- 1) .* (xi_squared * sss - dimension(q)) / 2
    a = 0.0

    return ρ .* weight .* (
        1.0 .+
        sss * u_dot_xi .+
        (sss^2 / 2) * (u_dot_xi .* u_dot_xi) .+
        a .+
        - (sss / 2) * u_squared
    )
end

function equilibrium(q::Quadrature, ρ, velocity, T)
    f = fill(0.0, length(q.weights))
    equilibrium!(q, f, ρ, velocity, T)
    return f
end

function hermite_equilibrium(q::Quadrature, f)
    ρ = density(q, f)
    j = momentum(q, f)
    u = j ./ ρ
    T = temperature(q, f, ρ, u)

    ξ = q.abscissae
    u_dot_ξ = dot(u, ξ)

    return ρ + ρ * u' * u + ρ * T * u' * u
end

function hermite_first_nonequilibrium(q::Quadrature, f)
end
