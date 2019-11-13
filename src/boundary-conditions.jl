abstract type BoundaryCondition end

"""
Apply a group of boundary conditions after streaming
"""
function apply!(boundary_conditions::Vector{<:BoundaryCondition}, q::Quadrature, f_new, f_old; time = 0.0)
    for boundary_condition in boundary_conditions
        apply!(boundary_condition, q, f_new, f_old, time = time)
    end
end

function apply!(bc::BoundaryCondition, q::Quadrature, f_new, f_old; time = 0.0)
    apply!(bc::BoundaryCondition, q::Quadrature, f_new, f_old)
end

"""
Often we need special rules when a boundary is placed in a specific direction

Currently we only support a few 2D boundary conditions
In the future we could add additional directions to support specific corners
"""
abstract type Direction end
struct North <: Direction end
struct East <: Direction end
struct South <: Direction end
struct West <: Direction end

# Halfway Bounce Back rules
struct BounceBack{D <: Direction, Ints} <: BoundaryCondition
    direction::D
    xs::Ints
    ys::Ints
end

function apply!(bc::BounceBack{<:North}, q::Quadrature, f_new, f_old)
    nx, ny, nf = size(f_new)
    for f_idx = 1:size(f_new, 3)
        # Bounce back the unkown populations
        opposite_f_idx = opposite(q, f_idx)

        for y_idx = bc.ys
            if y_idx + q.abscissae[2, opposite_f_idx] <= ny
                continue
            end
            for x_idx = bc.xs
                f_new[x_idx, y_idx, f_idx] = f_old[x_idx, y_idx, opposite_f_idx]
            end
        end
    end
end

function apply!(bc::BounceBack{<:South}, q::Quadrature, f_new, f_old)
    nx, ny, nf = size(f_new)
    for f_idx = 1:size(f_new, 3)
        opposite_f_idx = opposite(q, f_idx)

        for y_idx = bc.ys
            if y_idx + q.abscissae[2, opposite_f_idx] > 0
                continue
            end
            for x_idx = bc.xs
                f_new[x_idx, y_idx, f_idx] = f_old[x_idx, y_idx, opposite_f_idx]
            end
        end
    end
end

function apply!(bc::BounceBack{<:East}, q::Quadrature, f_new, f_old)
    nx, ny, nf = size(f_new)
    for f_idx = 1:size(f_new, 3)
        # Bounce back the unkown populations
        opposite_f_idx = opposite(q, f_idx)

        for x_idx = bc.xs
            if x_idx + q.abscissae[1, opposite_f_idx] <= nx
                continue
            end
            for y_idx = bc.ys
                f_new[x_idx, y_idx, f_idx] = f_old[x_idx, y_idx, opposite_f_idx]
            end
        end
    end
end
function apply!(bc::BounceBack{<:West}, q::Quadrature, f_new, f_old)
    nx, ny, nf = size(f_new)
    for f_idx = 1:size(f_new, 3)
        # Bounce back the unkown populations
        opposite_f_idx = opposite(q, f_idx)

        for x_idx = bc.xs
            if x_idx + q.abscissae[1, opposite_f_idx] > 0
                continue
            end
            for y_idx = bc.ys
                f_new[x_idx, y_idx, f_idx] = f_old[x_idx, y_idx, opposite_f_idx]
            end
        end
    end
end



"""
The MovingWall boundary condition applies a bounce back rule and adds a velocity
to each distribution that is bounced back
"""
struct MovingWall{D <: Direction, Ints} <: BoundaryCondition
    direction::D
    xs::Ints
    ys::Ints
    u::Vector{Float64}
    ρ::Float64
    T::Float64
end
MovingWall(direction, xs, ys, u) = MovingWall(direction, xs, ys, u, 1.0, 1.0)

function apply!(bc::MovingWall{<:North}, q::Quadrature, f_new, f_old)
    cs = q.speed_of_sound_squared

    a_1_eq = equilibrium_coefficient(Val{1}, q, bc.ρ, bc.u, bc.T)
    nx, ny, nf = size(f_new)
    for f_idx = 1:size(f_new, 3)
        H_i = hermite(Val{1}, q.abscissae[:, f_idx], q)
        a_1 = q.weights[f_idx] * cs * dot(a_1_eq, H_i)

        opposite_f_idx = opposite(q, f_idx)

        for y_idx = 1:ny
            if y_idx + q.abscissae[2, opposite_f_idx] <= ny
                continue
            end
            for x_idx = 1:nx
                f_new[x_idx, y_idx, f_idx] = f_old[x_idx, y_idx, opposite_f_idx] + 2 * a_1
            end
        end

    end
end
