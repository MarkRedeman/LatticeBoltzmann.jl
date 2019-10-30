
# Halfway bounce back rules

function bounce_back_from_bottom!(q::Quadrature, f_new, f_old, xs, ys)
    for f_idx = 1:size(f_new, 3)
        # Bounce back the unkown populations
        if q.abscissae[2, f_idx] > 0
            f_new[xs, ys, f_idx] = f_old[xs, ys, opposite(q, f_idx)]
        end
    end
end

function bounce_back_from_top!(q::Quadrature, f_new, f_old, xs, ys)
    for f_idx = 1:size(f_new, 3)
        # Bounce back the unkown populations
        if q.abscissae[2, f_idx] < 0
            f_new[xs, ys, f_idx] = f_old[xs, ys, opposite(q, f_idx)]
        end
    end
end

function bounce_back_from_left!(q::Quadrature, f_new, f_old, xs, ys)
    for f_idx = 1:size(f_new, 3)
        # Bounce back the unkown populations
        if q.abscissae[1, f_idx] > 0
            f_new[xs, ys, f_idx] = f_old[xs, ys, opposite(q, f_idx)]
        end
    end
end

function bounce_back_from_right!(q::Quadrature, f_new, f_old, xs, ys)
    for f_idx = 1:size(f_new, 3)
        # Bounce back the unkown populations
        if q.abscissae[1, f_idx] < 0
            f_new[xs, ys, f_idx] = f_old[xs, ys, opposite(q, f_idx)]
        end
    end
end

function bounce_back_from_top!(q::Quadrature, f_new, f_old, xs, ys, u::Vector{Float64})
    ρ_wall = 1.0
    for f_idx = 1:size(f_new, 3)
        # Bounce back the unkown populations
        if q.abscissae[2, f_idx] < 0
            # f_new[xs, ys, f_idx] = f_old[xs, ys, opposite(q, f_idx)] .+ 0.5 * (1 / q.speed_of_sound_squared) * q.abscissae[1, f_idx] * u[1]
            moving_wall = 2 * ρ_wall * q.speed_of_sound_squared * q.weights[f_idx] * q.abscissae[1, f_idx] * u[1]

            # @show moving_wall = 0.5 * (1 / q.speed_of_sound_squared) * q.abscissae[1, f_idx] * u[1]

            f_new[xs, ys, f_idx] = f_old[xs, ys, opposite(q, f_idx)] .+ moving_wall
        end
    end
end
