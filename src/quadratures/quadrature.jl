abstract type Quadrature end
abstract type Lattice end

function initialize(quadrature::Quadrature, lattice::Lattice)

end

density(q::Quadrature, f) = sum(f, dims=3)
function momentum(q::Quadrature, f)
    jx = sum([
        f[x, y, f_idx] * q.abscissae[1, f_idx] for x = 1 : size(f, 1), y = 1 : size(f, 2), f_idx = 1 : size(f, 3)
    ], dims = 3)
    jy = sum([
        f[x, y, f_idx] * q.abscissae[2, f_idx] for x = 1 : size(f, 1), y = 1 : size(f, 2), f_idx = 1 : size(f, 3)
    ], dims = 3)

    # jx = sum(f[:,:,[6, 7, 8]], dims=3) .- sum(f[:,:,[2, 3, 4]], dims=3);
    # jy = sum(f[:,:,[2, 8, 9]], dims=3) .- sum(f[:,:,[4, 5, 6]], dims=3);

    jx, jy
end

function total_energy(q::Quadrature, f)
    return sum([
        f[x, y, f_idx] * (q.abscissae[:, f_idx]' * q.abscissae[:, f_idx]) for x = 1 : size(f, 1), y = 1 : size(f, 2), f_idx = 1 : size(f, 3)
    ], dims = 3)
end

function internal_energy(q::Quadrature, f, u)

end

function kinetic_energy(q::Quadrature, f, u)

end

function equilibrium(quadrature::Quadrature, f_in)
    # Density
    f_ρ = density(quadrature, f_in)

    # Momentum
    jx, jy = momentum(quadrature, f_in)

    # Velocity componetns
    u = jx ./ f_ρ
    v = jy ./ f_ρ

    feq = equilibrium(quadrature, f_ρ, u, v, 1.0);

    # Compute equilibrium distribution
    return feq
end

function equilibrium(q::Quadrature, ρ, u, v, T)
    f = zeros(size(ρ,1), size(ρ,2), length(q.weights));

    u_squared = u.^2 + v.^2
    for idx = 1:length(q.weights)
        cs = q.abscissae[1, idx] .* u .+ q.abscissae[2, idx] .* v
        # cs = dot(abscissae[:, idx], u)


        # f[:, :, idx] = ρ .* weights[idx] .* (1.0 .+ 3.0 * cs .+ 4.5 * (cs .* cs) .- 1.5 * u_squared)
        f[:, :, idx] = _equilibrium(ρ, q.weights[idx], cs, u_squared)
    end

    return f
end

function equilibrium!(q::Quadrature, f, ρ, velocity, T)
    u = velocity[1]
    v = velocity[2]

    u_squared = u.^2 .+ v.^2

    for idx = 1:length(q.weights)
        cs = q.abscissae[1, idx] .* u .+ q.abscissae[2, idx] .* v


        f[idx] = _equilibrium(ρ, q.weights[idx], cs, u_squared)
        # ρ .* q.weights[idx] .* (1.0 .+ 3.0 * cs .+ 4.5 * (cs .* cs) .- 1.5 * u_squared)
    end

    jx = sum([f[f_idx] * q.abscissae[1, f_idx] for f_idx = size(f, 3)], dims = 3)
    jy = sum([f[f_idx] * q.abscissae[2, f_idx] for f_idx = size(f, 3)], dims = 3)

    return f
end


function _equilibrium(ρ, weight, cs, u_squared)
    return ρ .* weight .* (1.0 .+ 3.0 * cs .+ 4.5 * (cs .* cs) .- 1.5 * u_squared)
end

function equilibrium(q::Quadrature, ρ, velocity, T)
    f = fill(0.0, length(q.weights))
    equilibrium!(q, f, ρ, velocity, T)
    return f
end
