module Example
module TaylorGreenVortex

using Plots

struct TaylorGreenVortexExample end

abstract type Quadrature end
struct D2Q9 <: Quadrature end
tNS(::D2Q9) = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36]
cxNS(::D2Q9) = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1]
cyNS(::D2Q9) = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1]
oppNS(::D2Q9) = [  1,   4,  5,  2,  3,    8,   9,   6,   7]


const cx = [1 0 -1  0 1 -1 -1  1 0];
const cy = [0 1  0 -1 1  1 -1 -1 0];
const original_order = [1, 7, 9, 3, 5, 8, 2, 4, 6]
const weights   = [4/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9]'
const abscissae = [
    0  -1  -1  -1   0   1  1  1  0
    0   1   0  -1  -1  -1  0  1  1
]



function equilibrium(::D2Q9, ρ, velocity, T)
    u = velocity[1]
    v = velocity[2]

    subexp1 = 1.0 .- 1.5*(u.^2+v.^2);

    f = fill(0.0, 9)
    f[1] = (1/9) *  ρ .* (subexp1 + 3 .* u+9/2 .* u.^2);
    f[2] = (1/9) *  ρ .* (subexp1 + 3 .* v+9/2 .* v.^2);
    f[3] = (1/9) *  ρ .* (subexp1 - 3*u+9/2*u.^2);
    f[4] = (1/9) *  ρ .* (subexp1 - 3 .* v+9/2*v.^2);

    f[5] = (1/36) * ρ .* (subexp1 + 3*( u+v)+9/2*( u+v).^2);
    f[6] = (1/36) * ρ .* (subexp1 + 3*(-u+v)+9/2*(-u+v).^2);
    f[7] = (1/36) * ρ .* (subexp1 + 3*(-u-v)+9/2*(-u-v).^2);
    f[8] = (1/36) * ρ .* (subexp1 + 3*( u-v)+9/2*( u-v).^2);

    f[9] = (4/9) * ρ .* (subexp1);

    return f
end

function equilibrium!(::D2Q9, f, ρ, velocity, T)
    # fill!(f, 1.0)
    u = velocity[1]
    v = velocity[2]

    subexp1 = 1.0 .- 1.5*(u.^2+v.^2);

    f[1] = (1/9) *  ρ .* (subexp1 + 3*u+9/2*u.^2);
    f[2] = (1/9) *  ρ .* (subexp1 + 3 .* v+9/2*v.^2);
    f[3] = (1/9) *  ρ .* (subexp1 - 3*u+9/2*u.^2);
    f[4] = (1/9) *  ρ .* (subexp1 - 3*v+9/2*v.^2);

    f[5] = (1/36) * ρ .* (subexp1 + 3*( u+v)+9/2*( u+v).^2);
    f[6] = (1/36) * ρ .* (subexp1 + 3*(-u+v)+9/2*(-u+v).^2);
    f[7] = (1/36) * ρ .* (subexp1 + 3*(-u-v)+9/2*(-u-v).^2);
    f[8] = (1/36) * ρ .* (subexp1 + 3*( u-v)+9/2*( u-v).^2);

    f[9] = (4/9) * ρ .*(subexp1);

    return
end

function stream!(f, f_next)
    f_next[:,:,1] = circshift(f[:,:,1],[ 0  1]);
    f_next[:,:,2] = circshift(f[:,:,2],[ 1  0]);
    f_next[:,:,3] = circshift(f[:,:,3],[ 0 -1]);
    f_next[:,:,4] = circshift(f[:,:,4],[-1  0]);
   
    f_next[:,:,5] = circshift(f[:,:,5],[ 1  1]);
    f_next[:,:,6] = circshift(f[:,:,6],[ 1 -1]);
    f_next[:,:,7] = circshift(f[:,:,7],[-1 -1]);
    f_next[:,:,8] = circshift(f[:,:,8],[-1  1]);

    f_next[:,:,9] = f[:,:,9];
end

"""
This stream function applies periodic streaming, that is distribution
functions f_i that stream out of the domain will be placed on the opposite
side of the domain
"""
function stream(f, f_new = copy(f))

    # f_new[:,:,1] = circshift(f[:,:,1],[ 0  1]);
    # f_new[:,:,2] = circshift(f[:,:,2],[ 1  0]);
    # f_new[:,:,3] = circshift(f[:,:,3],[ 0 -1]);
    # f_new[:,:,4] = circshift(f[:,:,4],[-1  0]);

    # f_new[:,:,5] = circshift(f[:,:,5],[ 1  1]);
    # f_new[:,:,6] = circshift(f[:,:,6],[ 1 -1]);
    # f_new[:,:,7] = circshift(f[:,:,7],[-1 -1]);
    # f_new[:,:,8] = circshift(f[:,:,8],[-1  1]);

    # f_new[:,:,9] = f[:,:,9];

    # return f_new


    lx, ly, lq = size(f)

    @inbounds for x = 1:lx, y = 1:ly, f_idx = 1:lq
        next_x, next_y = stream_periodically_to(x, y, lx, ly, f_idx)

        f_new[next_x, next_y, f_idx] = f[x, y, f_idx]
    end

    return f_new
end

function equilibrium(rho, u, u_squared, idx::Int)
    cs = dot(abscissae[:, idx], u)

    return rho * weights[idx] .* (1.0 + 3.0 * cs + 4.5 * (cs .* cs) - 1.5 * u_squared)
end

function equilibrium(ρ,u,v)
    f = zeros(size(ρ,1), size(ρ,2), 9);

    subexp1 = 1 .- 1.5*(u.^2+v.^2);

    f[:, :, 1] = (1/9) .*  ρ .* (subexp1 + 3 .* u .+ 9 ./ 2 .* u.^2);
    f[:, :, 2] = (1/9) .*  ρ .* (subexp1 + 3 .* v .+ 9 ./ 2 .* v.^2);
    f[:, :, 3] = (1/9) .*  ρ .* (subexp1 - 3 .* u .+ 9 ./ 2 .* u.^2);
    f[:, :, 4] = (1/9) .*  ρ .* (subexp1 - 3 .* v .+ 9 ./ 2 .* v.^2);

    f[:, :, 5] = (1/36) .* ρ .* (subexp1 + 3 .* ( u .+ v) .+ 9 ./ 2 .* ( u .+ v).^2);
    f[:, :, 6] = (1/36) .* ρ .* (subexp1 + 3 .* (-u .+ v) .+ 9 ./ 2 .* (-u .+ v).^2);
    f[:, :, 7] = (1/36) .* ρ .* (subexp1 + 3 .* (-u .- v) .+ 9 ./ 2 .* (-u .- v).^2);
    f[:, :, 8] = (1/36) .* ρ .* (subexp1 + 3 .* ( u .- v) .+ 9 ./ 2 .* ( u .- v).^2);

    f[:, :, 9] = (4/9) .* ρ .* (subexp1);

    u_squared = u.^2 + v.^2
    for idx = 1:9
        cs = abscissae[1, idx] .* u .+ abscissae[2, idx] .* v


        f[:, :, idx] = ρ .* weights[idx] .* (1.0 .+ 3.0 * cs .+ 4.5 * (cs .* cs) .- 1.5 * u_squared)
    end

    return f
end

function collide(f_in, τ)
    f_out = copy(f_in)

    # Density
    f_ρ = sum(f_in, dims = 3)

    # Momentum components
    jx = sum(f_in[:,:,[1, 5, 8]], dims=3) .- sum(f_in[:,:,[3, 6, 7]], dims=3);
    jy = sum(f_in[:,:,[2, 5, 6]], dims=3) .- sum(f_in[:,:,[4, 7, 8]], dims=3);

    # Velocity componetns
    u = jx ./ f_ρ
    v = jy ./ f_ρ

    # @show u.^2 .+ v.^2
    # @show u

    # Compute equilibrium distribution
    feq = equilibrium(f_ρ, u, v);

    @show sum(f_in), sum(feq), sum(jx.^2 .+ jy.^2)

    return feq
    # f_out = feq
    f_out = (1 - 1 / τ) * f_in + (1 / τ) * feq;
    return f_out
    return f_in
    # Collision step

    return f_out
end

"""
Choose the next indices which should be streamed to depending on the given
 x and y index and the direction index.
We use the global abscissae variable to determine the direction and make
 sure that the indices are bounded
"""
function stream_periodically_to(x, y, lx, ly, f_idx)
    # Note: to do circshift: we have to subtract
    next_x = x - cx[f_idx]
    if next_x > lx
        next_x -= lx
    elseif next_x < 1
        next_x += lx
    end

    next_y = y - cy[f_idx]
    if next_y > ly
        next_y -= ly
    elseif next_y < 1
        next_y += ly
    end

    return next_x, next_y
end

density(x, y) = 1.0
velocity(x, y) = [
    cos(x) * sin(y),
    sin(x) * cos(y)
]

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
function process(f_in, t, stats)
    if (mod(t, 5) == 0)
        f_ρ = sum(f_in, dims = 3)

        # Momentum components
        jx = sum(f_in[:,:,[1, 5, 8]], dims=3) .- sum(f_in[:,:,[3, 6, 7]], dims=3);
        jy = sum(f_in[:,:,[2, 5, 6]], dims=3) .- sum(f_in[:,:,[4, 7, 8]], dims=3);

        # Velocity componetns
        u = jx ./ f_ρ
        v = jy ./ f_ρ

        # @show u.^2 .+ v.^2
        # @show u
        plot(
            contour(f_ρ[:, :, 1], size=(600, 600), fill=true),
            contour(u[:, :, 1].^2 .+ v[:, :, 1].^2, size=(600, 600), fill=true),
            size=(1200, 600)
        )
        gui()
    end
end

function siumlate(::TaylorGreenVortexExample;)
    τ = 1.0;

    N = 2^6 + 1
    # N = 129
    grid_size = (N, N)
    q = 9

    # initialize
    f_out = fill(0.0, grid_size..., 9)
    f_in = copy(f_out)

    f_ρ = fill(1.0, grid_size..., 1)
    f_u = fill(0.0, grid_size..., 2)
    f_T = fill(0.0, grid_size..., 1)

    # initial fields
    grid = [(x = x, y = y) for x in range(0, 2pi, length = N), y in range(0, 2pi, length = N)]

    density_field = [
        density(x, y) for x in range(0, 2pi, length = N), y in range(0, 2pi, length = N)
    ]

    velocity_field = [
        # velocity(x, y) for x in range(0, 2pi, length = N), y in range(0, 2pi, length = N)
        # [0.01 0.0] for x in range(0, 2pi, length = N), y in range(0, 2pi, length = N)
        0.005 * rand(2) for x in range(0, 2pi, length = N), y in range(0, 2pi, length = N)
    ]

    # set initial conditions
    feq = zeros(9)
    for x_idx = 1 : grid_size[1], y_idx = 1 : grid_size[2]
        # @views f = f_in[x_idx, y_idx, :]
        ρ = density_field[x_idx, y_idx]
        u = velocity_field[x_idx, y_idx]
        T = 1.0

        # f_in[x_idx, y_idx, :] = equilibrium(D2Q9(), ρ, u, T)
        f_out[x_idx, y_idx, :] = equilibrium(D2Q9(), ρ, u, T)


        f = f_out[x_idx, y_idx, :]
        ρ = sum(f)
        u = zeros(2)
        u[1] = (sum(f[[1 5 8]]) - sum(f[[3 6 7]])) ./ ρ
        u[2] = (sum(f[[2 5 6]]) - sum(f[[4 7 8]])) ./ ρ
        T = 1.0

        f_ρ[x_idx, y_idx] = ρ
        f_u[x_idx, y_idx, :] = u
        f_T[x_idx, y_idx] = T
    end

    contour(f_u[:, :, 1].^2 + f_u[:, :, 2].^2, size=(600, 600), fill=true)
    gui()

    stats = []

    @inbounds for t = 1 : 4N
        f_in = stream(f_out)

        process!(f_in, t, stats)

        # collide
        # collide(lattice, collision_model)
        f_out = collide(f_in, τ)


        continue;

        for x_idx = 1 : grid_size[1], y_idx = 1 : grid_size[2]
            # Get views (for performance?)
            f = f_in[x_idx, y_idx, :]

            # Compute moments
            ρ = sum(f)
            u[1] = (sum(f[[1 5 8]]) - sum(f[[3 6 7]])) ./ ρ
            u[2] = (sum(f[[2 5 6]]) - sum(f[[4 7 8]])) ./ ρ
            T = 1.0

            f_ρ[x_idx, y_idx] = ρ
            f_u[x_idx, y_idx, :] = u
            f_T[x_idx, y_idx] = T

            # f_in[x_idx, y_idx, :] = equilibrium(D2Q9(), ρ, u, T)

            feq = equilibrium(D2Q9(), ρ, u, T)
            # equilibrium!(D2Q9(), feq, ρ, u, T)
            #
            if (x_idx == 3 && y_idx == 3)
                @show f_in[x_idx, y_idx, 2], f_out[x_idx, y_idx, 2], f[2], feq[2]
                # @show ρ, u
                # @show (1 - 1 / τ) * f[2] + (1 / τ) * feq[2];
                # @show f[2]
                # @show feq[2];
            end
            f_out[x_idx, y_idx, :] = feq
            # f_out[x_idx, y_idx, :] = @. (1 - 1 / τ) * f + (1 / τ) * feq;
        end

        # stream
        # stream(lattice)
        f_in = stream(f_out)

        # post processing
        # @show sum(f_in), sum(f_out)
        if (mod(t, 1) == 0)
            plot(
                contour(f_u[:, :, 1].^2 + f_u[:, :, 2].^2, size=(600, 600), fill=true),
                # contour(f_u[:, :, 1], size=(600, 600)),
                # contour(f_u[:, :, 2], size=(600, 600)),
                # size=(3 * 600, 600),
                # layout=(1, 3)
            )
            gui()
        end
    end


    f_in, f_ρ, f_u, f_T
end

example = TaylorGreenVortexExample()

@time result = siumlate(example)

end
end
