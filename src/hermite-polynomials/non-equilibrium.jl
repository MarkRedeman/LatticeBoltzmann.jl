
function non_equilibrium(
    q::Quadrature,
    ρ::Float64,
    u::Vector{Float64},
    T::Float64
)
    fs = zeros(size(ρ,1), size(ρ,2), length(q.weights));
    f = zeros(length(q.weights));
    for x_idx = 1 : size(ρ, 1), y_idx = 1 : size(ρ, 2)
        equilibrium!(q, ρ[x_idx, y_idx], u[x_idx, y_idx, :], T[x_idx, y_idx], f)
        fs[x_idx, y_idx, :] .= f
    end
    return fs
end
