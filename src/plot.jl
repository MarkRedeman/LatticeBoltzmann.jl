function plot_density(fIn::Array{Float64, 3})
    lx = size(fIn, 2)
    ly = size(fIn, 3)

    # MACROSCOPIC VARIABLES
    rho = zeros(lx, ly)
    for x = 1:lx, y = 1:ly
        rho[x, y] = sum(fIn[:, x, y])
    end

    contour(
        1:lx,
        1:ly,
        rho,
        fill = true
    )
    # plot(layer(x=1:lx, y=1:ly, z=rho, Geom.contour))
end

function plot_speed(fIn::Array{Float64, 3})
    lx = size(fIn, 2)
    ly = size(fIn, 3)

    rho = zeros(lx, ly)
    ux = zeros(lx, ly)
    uy = zeros(lx, ly)

    # MACROSCOPIC VARIABLES
    for x = 1:lx, y = 1:ly
        rho[x, y] = sum(fIn[:, x, y])
        ux[x, y] = sum(fIn[:, x, y] .* cxNS) / rho[x, y]
        uy[x, y] = sum(fIn[:, x, y] .* cyNS) / rho[x, y]
    end

    contour(
        1:lx,
        1:ly,
        ux.^2 + uy.^2,
        fill = true
    )
end
