struct MRT{
    Force,
    VT <: AbstractVector{Float64},
    HST <: AbstractVector{Array{T,1} where T},
    AST <: AbstractVector{Array{T,1} where T}
} <: CollisionModel
    τs::VT

    # We will be using Hermite polynomials to compute the corresponding coefficients
    H0::Float64
    Hs::HST
    # Keep higher order coefficients allocated
    As::AST

    force::Force
end

function MRT(q::Quadrature, τ::Float64, force = nothing)
    N = round(Int, order(q) / 2)
    MRT(q, fill(τ, N))
end

function MRT(q::Quadrature, τ_s::Float64, τ_a::Float64, force = nothing)
    N = round(Int, order(q) / 2)
    MRT(q, repeat([τ_s, τ_a], outer = N)[1:N])
end

function MRT(q::Quadrature, τs::VT, force = nothing) where { VT <: AbstractVector{Float64} }
    N = round(Int, order(q) / 2)
    Hs = [[hermite(Val{n}, q.abscissae[:, i], q) for i = 1:length(q.weights)] for n = 1:N]

    MRT(τs, 1.0, Hs, copy(Hs), force)
end

function CollisionModel(cm::Type{<:MRT}, q::Quadrature, problem::FluidFlowProblem)
    τ = q.speed_of_sound_squared * lattice_viscosity(problem) + 0.5
    τs = fill(τ, order(q))

    if has_external_force(problem)
        force_field = (x_idx, y_idx, t) -> lattice_force(problem, x_idx, y_idx, t)

        return MRT(q, τs, force_field)
    end

    return MRT(q, τs)
end

function collide!(collision_model::MRT{Force}, q, f_in, f_out; time = 0.0) where {Force}
    collide_mrt!(collision_model, q, f_in, f_out, time = time)
end

# NOTE: we can do some clever optimizations where we check the value of τ_n
# If it is equal to 1, then we don't need to compute the nth hermite coefficient
# of f
function collide_mrt!(collision_model::MRT{Force}, q, f_in, f_out; time = 0.0) where {Force}
    cs = q.speed_of_sound_squared

    τs = collision_model.τs

    D = dimension(q)
    N = div(lbm.order(q), 2)

    # Compute (get!?) the hermite polynomials for this quadrature
    Hs = collision_model.Hs

    nx, ny, nf = size(f_in)

    if !(Force <: Nothing)
        F = zeros(dimension(q))
    end

    f = Array{Float64}(undef, nf)
    u = zeros(dimension(q))
    @inbounds for x_idx = 1:nx, y_idx = 1:ny
        @inbounds for f_idx = 1:nf
            f[f_idx] = f_in[x_idx, y_idx, f_idx]
        end

        ρ = density(q, f)
        velocity!(q, f, ρ, u)
        # T = temperature(q, f, ρ, a_f[1], a_f[2])
        T = 1.0

        if !(Force <: Nothing)
            F .= collision_model.force(x_idx, y_idx, time)

            u += τs[2] * F
        end

        # NOTE: we don't need to compute the 0th and 1st coefficient as these are equal
        # to a_f[0] and a_f[1]
        a_eq = [equilibrium_coefficient(Val{n}, q, ρ, u, T) for n = 1:N]

        # NOTE: we could optimize this by only computing upto n = 2 when τ = 1.0
        a_f = [sum([f[idx] * Hs[n][idx] for idx = 1:length(q.weights)]) for n = 1:N]

        a_coll = [(1 - 1 / τs[n]) .* a_f[n] .+ (1 / τs[n]) .* a_eq[n] for n = 1:N]

        @inbounds for f_idx = 1:nf
            f_out[x_idx, y_idx, f_idx] =
                q.weights[f_idx] * (
                    ρ +
                    cs * ρ * sum(u .* Hs[1][f_idx]) +
                    sum([cs^n * dot(a_coll[n], Hs[n][f_idx]) / (factorial(n)) for n = 2:N])
                )
        end
    end
    return
end
