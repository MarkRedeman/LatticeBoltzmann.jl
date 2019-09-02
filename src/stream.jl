# const cx = [1 0 -1  0 1 -1 -1  1 0];
# const cy = [0 1  0 -1 1  1 -1 -1 0];
const cx = [0  -1  -1  -1   0   1  1  1  0]
const cy = [0   1   0  -1  -1  -1  0  1  1]
"""
This stream function applies periodic streaming, that is distribution
functions f_i that stream out of the domain will be placed on the opposite
side of the domain
"""
function stream(f, f_new = copy(f))
    lx, ly, lq = size(f)

    @inbounds for x = 1:lx, y = 1:ly, f_idx = 1:lq
        next_x, next_y = stream_periodically_to(x, y, lx, ly, f_idx)

        f_new[next_x, next_y, f_idx] = f[x, y, f_idx]
    end

    return f_new
end

function stream(quadrature, f, f_new = copy(f))
    lx, ly, lq = size(f)

    @inbounds for x = 1:lx, y = 1:ly, f_idx = 1:lq
        next_x, next_y = stream_periodically_to(quadrature, x, y, lx, ly, f_idx)

        f_new[next_x, next_y, f_idx] = f[x, y, f_idx]
    end

    return f_new
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

function stream_periodically_to(q::D2Q9, x, y, lx, ly, f_idx)
    # Note: to do circshift: we have to subtract
    next_x = x - abscissae[1, f_idx]
    if next_x > lx
        next_x -= lx
    elseif next_x < 1
        next_x += lx
    end

    next_y = y - abscissae[2, f_idx]
    if next_y > ly
        next_y -= ly
    elseif next_y < 1
        next_y += ly
    end

    return next_x, next_y
end
