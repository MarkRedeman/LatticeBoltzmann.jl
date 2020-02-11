# Use F.2 to generate them
function hermite_based_equilibrium(q, f)
    ρ = sum(f)
    u = velocity(q, f)
    T = temperature(q, f, ρ, u)

    hermite_based_equilibrium(q, ρ, u, T)
end

function hermite_based_equilibrium!(q::Q, ρ, u, T, f) where {Q <: Quadrature}
    N = div(LatticeBoltzmann.order(q), 2)

    # NOTE: we skip H_0 since its one (also because of Julia's 1 based index....)
    #Hs = [[hermite(n, q.abscissae[:, i], q) for i = 1:length(q.weights)] for n = 1:N]

    cs = 1 / q.speed_of_sound_squared
    D = dimension(q)
    δ(α, β) = α == β ? 1 : 0

    a0 = ρ
    H1 = [hermite(Val{1}, q.abscissae[:, i], q) for i in 1:length(q.weights)]
    a1 = u

    Hs = [[hermite(Val{n}, q.abscissae[:, i], q) for i in 1:length(q.weights)] for n in 1:N]
    a_eq = [equilibrium_coefficient(Val{n}, q, ρ, u, T) for n in 1:N]

    @inbounds for f_idx in 1:length(f)
        f[f_idx] =
            q.weights[f_idx] *
            (ρ + sum([dot(a_eq[n], Hs[n][f_idx]) / (factorial(n) * cs^n) for n in 1:N]))
    end
    return
end

# equilibrium(Val{1}, ρ, u, T)

function equilibrium_coefficient(::Type{Val{0}}, q::Quadrature, ρ, u, T)
    return ρ
end

function equilibrium_coefficient(::Type{Val{1}}, q::Quadrature, ρ, u, T)
    return ρ * u
end

function equilibrium_coefficient(::Type{Val{2}}, q::Quadrature, ρ, u, T)
    cs = 1 / q.speed_of_sound_squared
    return ρ * (u * u' + cs * (T - 1) * I)
end
function equilibrium_coefficient(::Type{Val{3}}, q::Quadrature, ρ, u, T)
    cs = 1 / q.speed_of_sound_squared
    D = dimension(q)
    return ρ * [
        u[a] * u[b] * u[c] + cs * (T - 1) * (
            u[a] * δ(b, c) + u[b] * δ(a, c) + u[c] * δ(a, b)
        )
        for a in 1:D, b in 1:D, c in 1:D
    ]
end
function equilibrium_coefficient(::Type{Val{4}}, q::Quadrature, ρ, u, T)
    cs = 1 / q.speed_of_sound_squared
    D = dimension(q)
    return ρ * [
        u[a] * u[b] * u[c] * u[d] +
        cs *
        (T - 1) *
        (
            u[a] * u[b] * δ(c, d) +
            u[a] * u[c] * δ(b, d) +
            u[a] * u[d] * δ(b, d) +
            u[b] * u[c] * δ(a, d) +
            u[b] * u[d] * δ(a, d) +
            u[c] * u[d] * δ(a, b)
        ) +
        cs^2 * (T - 1)^2 * (δ(a, b) * δ(c, d) + δ(a, c) * δ(b, d) + δ(a, d) * δ(b, c))
        for a in 1:D, b in 1:D, c in 1:D, d in 1:D
    ]
end

"""
Compute the temperature from hermite coefficients
"""
function temperature(
    q::Q,
    f::VT,
    a_0::VT,
    a_1::VT,
    a_2::MT,
) where {Q <: Quadrature, T <: Real, VT <: AbstractVector{T}, MT <: AbstractMatrix{T}}
    D = dimension(q)
    P = Array{T}(undef, D, D)
    ρ = a_0
    u = a_1 / ρ
    P = q.speed_of_sound_squared * a_2 - ρ * (u * u' - I(2))

    return tr(P) / (D * ρ)
end
