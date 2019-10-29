
# Halfway bounce back rules

function bounce_back_from_bottom!(q::Quadrature, f_new, f_old, xs, ys)
    for f_idx = 1:size(f_new, 3)
        if q.abscissae[2, f_idx] > 0
            f_new[xs, ys, f_idx] = f_old[xs, ys, opposite(q, f_idx)]
        end
    end
end

function bounce_back_from_top!(q::Quadrature, f_new, f_old, xs, ys)
    for f_idx = 1:size(f_new, 3)
        if q.abscissae[2, f_idx] < 0
            f_new[xs, ys, f_idx] = f_old[xs, ys, opposite(q, f_idx)]
        end
    end
end
