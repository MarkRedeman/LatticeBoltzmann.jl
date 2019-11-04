"""
This stream function applies periodic streaming, that is distribution
functions f_i that stream out of the domain will be placed on the opposite
side of the domain
"""
function stream(quadrature::Quadrature, f, f_new = copy(f))
    lx, ly, lq = size(f)

    @inbounds for x = 1:lx, y = 1:ly, f_idx = 1:lq
        next_x, next_y = stream_periodically_from(quadrature, x, y, lx, ly, f_idx)

        f_new[next_x, next_y, f_idx] = f[x, y, f_idx]
    end

    return f_new
end

stream!(q::Quadrature; f_new, f_old) = stream!(q, f_old, f_new)
function stream!(quadrature::Quadrature, f, f_new)
    lx, ly, lq = size(f)

    @inbounds for x = 1:lx, y = 1:ly, f_idx = 1:lq
        # next_x, next_y = stream_periodically_from(quadrature, x, y, lx, ly, f_idx)

        # f_new[next_x, next_y, f_idx] = f[x, y, f_idx]

        # Gather
        from_x, from_y = stream_periodically_to(quadrature, x, y, lx, ly, f_idx)

        f_new[x, y, f_idx] = f[from_x, from_y, f_idx]
    end

    return
end

"""
Choose the next indices which should be streamed to depending on the given
 x and y index and the direction index.
We use the global abscissae variable to determine the direction and make
 sure that the indices are bounded
"""
function stream_periodically_from(q::Quadrature, x, y, lx, ly, f_idx)
    # Note: to do circshift: we have to subtract
    next_x = x + q.abscissae[1, f_idx]
    if next_x > lx
        next_x -= lx
    elseif next_x < 1
        next_x += lx
    end

    next_y = y + q.abscissae[2, f_idx]
    if next_y > ly
        next_y -= ly
    elseif next_y < 1
        next_y += ly
    end

    return next_x, next_y
end


"""
Choose the next indices which should be streamed to depending on the given
 x and y index and the direction index.
We use the abscissae variable to determine the direction and make
 sure that the indices are bounded
"""
function stream_periodically_to(q::Quadrature, x, y, lx, ly, f_idx)
    next_x = x - q.abscissae[1, f_idx]
    next_y = y - q.abscissae[2, f_idx]

    return mod1(next_x, lx), mod1(next_y, ly)
end
