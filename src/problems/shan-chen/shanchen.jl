#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# shanChen.m: Multi-component fluid, using a LB method,
#   based on the Shan-Chen model
# [X.Shan and H.Chen, http://dx.doi.org/10.1103/PhysRevE.47.1815].
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Lattice Boltzmann sample, written in Matlab
# Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani, Jonas Latt
# Address: EPFL-STI-LIN Station 9
# E-mail: orestis.malaspinas@epfl.ch
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module ShanChen

# using Gadfly
using Plots

# GENERAL FLOW CONSTANTS
const G = -1.2;  # Amplitude of the molecular interaction force

const omega1 = 1.;  # Relaxation parameter for fluid 1
const omega2 = 1.;  # Relaxation parameter for fluid 2

const tPlot  = 100;  # iterations between successive graphical outputs

# D2Q9 LATTICE CONSTANTS
const tNS   = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
const cxNS  = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1];
const cyNS  = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1];
const oppNS = [  1,   4,  5,  2,  3,    8,   9,   6,   7];

dot(a, b) = sum(a .* b);

function run(lx::Int = 128, ly::Int = 128, maxT::Int = 4000)
    # GENERAL FLOW CONSTANTS
    drho = 0.001;
    delta_rho = -drho * (1.0 .- 2.0 * rand(lx, ly));

    # INITIAL CONDITION FOR BOTH DISTRIBUTION FUNCTIONS: (T=0) ==> TIn(i) = t(i)
    fIn = zeros(9, lx, ly)
    gIn = zeros(9, lx, ly)
    for i = 1:9, x = 1:lx, y = 1:ly
        fIn[i, x, y] = tNS[i] .* (1.0 .+ delta_rho[x, y]);
        gIn[i, x, y] = tNS[i] .* (1.0 .- delta_rho[x, y]);
    end

    # return fIn, gIn
    # STREAMING STEP
    # Stream fOut to new distributions
    return simulate(fIn, gIn, maxT)
end

function simulate(fIn::Array{Float64, 3}, gIn::Array{Float64, 3}, maxT::Int64 = 4000)
    lx = size(fIn, 2)
    ly = size(fIn, 3)

    rho1 = sum(fIn, dims=1)[1, :, :]
    rho2 = sum(gIn, dims=1)[1, :, :]

    # MAIN LOOP (TIME CYCLES)
    Gomega1 = G / omega1;
    Gomega2 = G / omega2;
    fOut = zeros(9, lx, ly)
    gOut = zeros(9, lx, ly)

    jx1 = zeros(1, lx, ly)
    jy1 = zeros(1, lx, ly)
    jx2 = zeros(1, lx, ly)
    jy2 = zeros(1, lx, ly)

    uTotX = zeros(lx, ly)
    uTotY = zeros(lx, ly)

    ρ = zeros(lx, ly)

    uTotX1 = zeros(lx, ly); #POTENTIAL CONTRIBUTION OF FLUID 2 ON 1
    uTotY1 = zeros(lx, ly);
    uTotX2 = zeros(lx, ly); #POTENTIAL CONTRIBUTION OF FLUID 2 ON 1
    uTotY2 = zeros(lx, ly);

    for cycle = 1:maxT
        rhoContrib1x = zeros(lx, ly);
        rhoContrib2x = zeros(lx, ly);

        rhoContrib1y = zeros(lx, ly);
        rhoContrib2y = zeros(lx, ly);

    # rho1 = sum(fIn, 1)[1, :, :]
    # rho2 = sum(gIn, 1)[1, :, :]

        # u_σ : 1×64×64 Array{Float64,3}:
        for x = 1 : lx, y = 1 : ly
            # MACROSCOPIC VARIABLES
            rho1[x, y] = sum(fIn[:, x, y])
            rho2[x, y] = sum(gIn[:, x, y])

            # Momentum
            jx1[1, x, y] = dot(cxNS, fIn[:, x, y])
            jy1[1, x, y] = dot(cyNS, fIn[:, x, y])

            jx2[1, x, y] = dot(cxNS, gIn[:, x, y])
            jy2[1, x, y] = dot(cyNS, gIn[:, x, y])

            # Density

            # (∑ ρ_σ / τ_σ)
            ρ[x, y] = rho1[x, y] * omega1 + rho2[x, y] * omega2;

            # u' = (∑ (ρ_σ u_σ)/ τ_σ) / (∑ ρ_σ / τ_σ)
            uTotX[x, y] = (jx1[1, x, y] * omega1 + jx2[1, x, y] * omega2) ./ ρ[x, y];
            uTotY[x, y] = (jy1[1, x, y] * omega1 + jy2[1, x, y] * omega2) ./ ρ[x, y];
            for i=2:9
                x_next = x - cxNS[i]
                y_next = y - cyNS[i]

                if x_next > lx
                    x_next -= lx
                elseif x_next < 1
                    x_next += lx
                end

                if y_next > ly
                    y_next -= ly
                elseif y_next < 1
                    y_next += ly
                end


                # Note!
                # Here we are computing from a non local density
                rhoContrib1x[x, y] = rhoContrib1x[x, y] + rho1[x_next, y_next] * tNS[i] * cxNS[i];
                rhoContrib1y[x, y] = rhoContrib1y[x, y] + rho1[x_next, y_next] * tNS[i] * cyNS[i];

                rhoContrib2x[x, y] = rhoContrib2x[x, y] + rho2[x_next, y_next] * tNS[i] * cxNS[i];
                rhoContrib2y[x, y] = rhoContrib2y[x, y] + rho2[x_next, y_next] * tNS[i] * cyNS[i];
            end

            uTotX1[x, y] = uTotX[x, y] - Gomega1 .* rhoContrib2x[x, y]; #POTENTIAL CONTRIBUTION OF FLUID 2 ON 1
            uTotY1[x, y] = uTotY[x, y] - Gomega1 .* rhoContrib2y[x, y];

            uTotX2[x, y] = uTotX[x, y] - Gomega2 .* rhoContrib1x[x, y]; #POTENTIAL CONTRIBUTION OF FLUID 2 ON 1
            uTotY2[x, y] = uTotY[x, y] - Gomega2 .* rhoContrib1y[x, y];

            # COLLISION STEP FLUID 1 AND 2
            # Compute the equilibrium distribution with the momentum change due to
            # the inter-molecular forces
            for i=1:9
                cuNS1 = (cxNS[i] * uTotX1[x, y] + cyNS[i] * uTotY1[x, y]);
                cuNS2 = (cxNS[i] * uTotX2[x, y] + cyNS[i] * uTotY2[x, y]);

                fEq   = rho1[x, y] .* tNS[i] .* ( 1 + 3.0 * cuNS1 + 4.5 * (cuNS1 .* cuNS1) - 1.5 * (uTotX1[x, y].^2 + uTotY1[x, y].^2));
                gEq   = rho2[x, y] .* tNS[i] .* ( 1 + 3.0 * cuNS2 + 4.5 * (cuNS2 .* cuNS2) - 1.5 * (uTotX2[x, y].^2 + uTotY2[x, y].^2));

                fOut[i,x,y]  = fIn[i,x,y] - omega1 .* (fIn[i,x,y] - fEq);
                gOut[i,x,y]  = gIn[i,x,y] - omega2 .* (gIn[i,x,y] - gEq);
            end
        end


        # STREAMING STEP FLUID 1 AND 2
        for i=1:9
            fIn[i, :, :] = stream(fOut[i, :, :], cxNS[i], cyNS[i])
            gIn[i, :, :] = stream(gOut[i, :, :], cxNS[i], cyNS[i])
        end

        # VISUALIZATION
        if (cycle % tPlot == 0)
            @show cycle
            plot_density(fIn)
            gui()
        end
    end
    return fIn, gIn
end


"""
    Periodic streaming
"""
function stream(fOut, cx, cy)
    circshift(fOut, [cx, cy]);
end

function equilibrium(i, ρ, u_x, u_y)
    cuNS1 = 3 * (cxNS[i] * u_x + cyNS[i] * u_y);

    return ρ .* tNS[i] .* (1 + cuNS1 + 1/2 * (cuNS1 .* cuNS1) - 3/2 * (u_x.^2 + u_y.^2));
end

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

# function plot_speed(ux::Array{Float64, 2}, uy::Array{Float64, 2}, rho::Array{Float64, 2})
#     lx = size(ux, 1)
#     ly = size(ux, 2)

#     u = reshape(sqrt(ux.^2 + uy.^2),lx,ly);
#     plot(layer(x=1:lx, y=1:ly, z=u, Geom.contour))
# end

end
