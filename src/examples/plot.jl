function streamline(j; step = round(Int, size(j, 1) / 5) )
    s = (1000, 500)
    velocity_field = contour(j[:, :, 1].^2 .+ j[:, :, 2].^2, cbar=true, fill=true, title="Momentum", size=s)
    N = size(j, 1)
    X = [i for i in range(1, size(j, 1), step = step), j in range(1, size(j, 2), step = step)]
    Y = [j for i in range(1, size(j, 1), step = step), j in range(1, size(j, 2), step = step)]
    # @show "process: ", u, v,

    quiver!(
        velocity_field,
        X, Y,
        quiver=(x, y) -> (j[x, y, 1] / sqrt(j[x, y, 1]^2 + j[x, y, 2]^2), j[x, y, 2] / sqrt(j[x, y, 1]^2 + j[x, y, 2]^2)) ,
        color="white",
    )

    return velocity_field
end
