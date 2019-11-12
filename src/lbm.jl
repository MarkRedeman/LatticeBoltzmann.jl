module lbm

using DataFrames
using BenchmarkTools
using TimerOutputs

include("quadratures/quadrature.jl")
include("stream.jl")
include("collision.jl")
include("problems/problems.jl")
include("process.jl")
include("stopping-criteria.jl")
include("boundary-conditions.jl")

export Quadrature, Lattice, D2Q4, D2Q5, D2Q9, D2Q17,
    opposite,
    initialize,
    density,
    momentum,
    pressure,
    total_energy,
    kinetic_energy,
    internal_energy,
    temperature,
    dimension,
    equilibrium,
    equilibrium!,
    hermite_equilibrium,
    hermite_first_nonequilibrium,
    hermite,
    stream,
    stream!,
    CollisionModel,
    SRT,
    SRT_Force,
    TRT,
    collide!,
    simulate

# https://github.com/MagB93/LatticeBoltzmann/blob/master/src/LatticeBoltzmann.jl
function siumlate(problem::InitialValueProblem, quadrature::Quadrature = D2Q9();
                  base = 200,
                  n_steps = base * problem.NX * problem.NX / (16 * 16),
                  should_process = true,
                  t_end = 1.0)
    # initialize
    f_out, collision_operator = initialize(quadrature, problem)
    f_in = copy(f_out)
    Δt = lbm.delta_t(problem)
    n_steps = round(Int, t_end / Δt)

    stats = process_stats()
    stop_criteria = StopCriteria(problem)

    @inbounds for t = 0:n_steps
        # if mod(t, round(Int, n_steps / 10)) == 0
        #     @show t, t / n_steps, t * Δt
        # end

        if (should_process && mod(t, 1) == 0)
            lbm.process!(
                problem,
                quadrature,
                f_in,
                t * Δt,
                stats,
                # should_visualize = (mod(t, 50) == 0)
                should_visualize = (mod(t, round(Int, n_steps / 20)) == 0)
                # should_visualize = true
            )
        end

        collide!(collision_operator, quadrature, time = t * Δt, problem = problem, f_old = f_in, f_new = f_out)

        stream!(quadrature, f_new = f_in, f_old = f_out)

        apply_boundary_conditions_after!(quadrature, problem, f_new = f_in, f_old = f_out, time = t * Δt)

        # check_stability(f_in) || return :unstable, f_in, stats       )
        if mod(t, 100) == 0
            if (should_stop!(stop_criteria, quadrature, f_in))
                @show "Stopping after $t steps out of $n_steps"
                lbm.process!(problem, quadrature, f_in, t * Δt, stats, should_visualize = should_process)
                return f_in, stats
            end
        end
    end

    lbm.process!(problem, quadrature, f_in, n_steps * Δt, stats, should_visualize = should_process)

    # @show stats

    f_in, stats
end

# export process!

# abstract type Model

# struct Problem{Model}
#     model::Model
#     # model::Model
#     # quadrature::Quadrature
#     # collision::Collision
#     N::Int64      # rename to points / cells / ...
#     N_iter::Int64 # rename to iterations?
# end

# # Put collision_model, relaxation_rate, simulate, initial_condition, here..

# # abstract Quadrature;
# # Distribution{D2Q9}
# # Distribution{D2Q17}
# # immutable Distribution{T, Quadrature}
# #     fs::Array{T, 1}
# # end
# # typealias Distribution{T} Distribution{Float64, Q}

# typealias Distributions Array{Float64, 3}
# typealias Distribution Array{Float64, 2}

# include("quadratures/D2Q9.jl")
# # include("quadratures/D2Q4.jl")

# density(f::Distribution) = sum(f)

# velocity(f::Distribution) = velocity(f, abscissae, original_order)
# function velocity{D}(f::Distribution, abscissae::Array{Int64, D}, order::Array{Int64, 1})
#     ρ, u = density_and_velocity(f, abscissae, order)
#     return u
# end

# density_and_velocity(f::Distribution) = density_and_velocity(f, abscissae, original_order)
# function density_and_velocity{D}(f::Distribution, abscissae::Array{Int64, D}, order::Array{Int64, 1})::Tuple{Float64, Array{Float64, 1}}
#     u = zeros(D)
#     ρ = 0.0

#     # Compute: ∑fᵢ and ∑fᵢξᵢ
#     for idx ∈ order
#         for d = 1:D
#             u[d] += abscissae[d, idx] * f[idx]
#         end
#         ρ += f[idx]
#     end

#     return ρ, u / ρ
# end

# equilibrium(f::Distribution) = equilibrium(density_and_velocity(f, abscissae, original_order)...)
# equilibrium{T}(ρ::T, u::Array{T, 1})::Distribution = equilibrium(ρ, u, dot(u, u))
# equilibrium{T}(ρ::T, u::Array{T, 1}, u_squared::T)::Distribution = [equilibrium(ρ, u, u_squared, idx) for idx = 1:9]'
# function equilibrium{T}(rho::T, u::Array{T, 1}, u_squared::T, idx::Int)::T
#     const cs = dot(abscissae[:, idx], u)

#     return rho * weights[idx] .* (1.0 + 3.0 * cs + 4.5 * (cs .* cs) - 1.5 * u_squared)
# end

# """
# By default we will be using a bgk collision
# """
# collide(f::Distribution, ω)::Distribution = bgk_collision(f, ω)
# collide(f::Distribution, ω, force)::Distribution = bgk_collision(f, ω, force)

# """
# Apply the most simple collision operator
# """
# function bgk_collision{T}(f::Distribution, ω::T)::Distribution
#     const ρ, u = density_and_velocity(f)

#     bgk_collision(f, ω, ρ, u)
# end

# """
# Apply the bgk collision operator after adding additional momentum
# to the velocity due to a force term
# """
# function bgk_collision{T}(f::Distribution, ω::T, force::Array{T, 1})::Distribution
#     const ρ, u = density_and_velocity(f)

#     bgk_collision(f, ω, ρ, u + force / ω)
# end

# function bgk_collision{T}(f::Distribution, ω::T, ρ::T, u::Array{T, 1})
#     const u_squared = dot(u, u)

#     for idx = 1:9
#         f[idx] = (1 - ω) * f[ if next_x > lx
#         next_x -= lx
#     elseif next_x < 1
#         next_x += lx
#     end

#     next_y = y - abscissae[2, f_idx]
#     if next_y > ly
#         next_y -= ly
#     elseif next_y < 1
#         next_y += ly
#     end

#     return next_x, next_y
# end



# """
#     step(lbm)
# Computes all individual steps needed in one iteration of the lattice boltzmann
# scheme.
#   1. collision operator ( only the BGK is implemented)
#   2. steaming of the distributions to neighbouring nodes
#   3. computation of the passed boundary conditions
#   4. computation of the macroskopic variables
#   5. equilibrium distribution function
# """
# function step!(grid::Grid, velset::Velocity_Set, collision::Collision,
#               stream::Array{Streaming, 1}, bound::Array{Boundary, 1})

#     compute_collision!(grid, collision)
#     compute_streaming!(grid, stream, velset)
#     compute_boundary!(grid, bound, velset)
#     compute_macro_var!(grid, velset)
#     compute_f_eq!(grid, velset)
# end

end
