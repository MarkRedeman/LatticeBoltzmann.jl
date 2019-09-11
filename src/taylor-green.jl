using Plots


module TaylorGreen

"""
Taylor Green Vortex in 2 dimensions
"""
abstract type TaylorGreenVortex{T} end

struct DecayingVortex{T<:Real} <: TaylorGreenVortex{T}
    a::T
    b::T
    A::T
    B::T
    Re::T # Reynolds number
    length::T
    speed::T

    function DecayingVortex(a, b, A, B, ν, length, speed)
        a * A + b * B != 0 ? throw(ArgumentError("The given vortex is not incompressible")) : new(a, b, A, B, ν, length, speed)
    end
end

struct StaticVortex{T<:Real} <: TaylorGreenVortex{T}
    a::T
    b::T
    A::T
    B::T
    Re::T # Reynolds number
    length::T
    speed::T

    function StaticVortex(a, b, A, B, ν, length, speed)
        a * A + b * B != 0 ? throw(ArgumentError("The given vortex is not incompressible")) : new(a, b, A, B, ν, length, speed)
    end
end

# Constructors
TaylorGreenVortex{T<:Real}(a::T, b::T, A::T, B::T, ν::T, length::T, speed::T) = DecayingVortex{T}(a, b, A, B, ν, length, speed)
TaylorGreenVortex{T<:Real}(a::T, b::T, A::T, B::T, ν::T) = DecayingVortex{T}(a, b, A, B, ν, 1., 1.)
TaylorGreenVortex() = TaylorGreenVortex{Float64}(1., 1., 1., -1., 1.)

velocity(t::TaylorGreenVortex, x, y) = velocity(t, x, y, 0.0)
function velocity(t::TaylorGreenVortex, x, y, time)
    return decay(t, time * t.length / t.speed) * (1 / t.speed) * [
        t.A * cos(t.a * t.length * x)sin(t.b * t.length * y),
        t.B * sin(t.a * t.length * x)cos(t.b * t.length * y)
    ]
end

decay(t::TaylorGreenVortex, time) = exp(- (t.a^2 + t.b^2) * (t.speed * t.length) / (t.Re) * time)
decay(t::StaticVortex, time) = 1
decay(t::DecayingVortex, time) = exp(- (t.a^2 + t.b^2) * (t.speed * t.length) / (t.Re) * time)

pressure(t::TaylorGreenVortex, x, y) = pressure(t, x, y, 0.0)
function pressure(t::TaylorGreenVortex, x, y, time)
    return (t.A^2 / (4 * t.speed^2)) * (cos(2 * t.a * t.length * x) + (t.a^2 / t.b^2) * cos(2 * t.b * t.length * y)) * decay(t, time * t.length / t.speed)
end

force(t::DecayingVortex, x, y) = 0
force(t::StaticVortex, x, y) = (t.a^2 + t.b^2) * (t.length^2 / t.Re) * velocity(t, x, y)


end


function taylorgreen(t, x, y, Nx, Ny, nu, rho0, u_max)
    kx = 2*pi/Nx;
    ky = 2*pi/Ny;
    td = 1/(nu*(kx*kx+ky*ky));

    u = -u_max*sqrt(ky/kx)*cos(kx*x).*sin(ky*y)*exp(-t/td);
    v =  u_max*sqrt(kx/ky)*sin(kx*x).*cos(ky*y)*exp(-t/td);
    P = -0.25*rho0*u_max*u_max*((ky/kx)*cos(2*kx*x) .+ (kx/ky)*cos(2*ky*y))*exp(-2*t/td);
    rho = rho0 .+ 3*P;

    return (rho, u, v, P)
end

function equilibrium!(feq, rho, u, v)
    # feq = zeros(size(rho,1),size(rho,2),9);

    subexp1 = 1.0 .- 1.5*(u.^2+v.^2);

    feq[:,:,1]=(1/9)*rho.*(subexp1 + 3*u+9/2*u.^2);
    feq[:,:,2]=(1/9)*rho.*(subexp1 + 3*v+9/2*v.^2);
    feq[:,:,3]=(1/9)*rho.*(subexp1 - 3*u+9/2*u.^2);
    feq[:,:,4]=(1/9)*rho.*(subexp1 - 3*v+9/2*v.^2);

    feq[:,:,5]=(1/36)*rho.*(subexp1 + 3*( u+v)+9/2*( u+v).^2);
    feq[:,:,6]=(1/36)*rho.*(subexp1 + 3*(-u+v)+9/2*(-u+v).^2);
    feq[:,:,7]=(1/36)*rho.*(subexp1 + 3*(-u-v)+9/2*(-u-v).^2);
    feq[:,:,8]=(1/36)*rho.*(subexp1 + 3*( u-v)+9/2*( u-v).^2);

    feq[:,:,9]=(4/9)*rho.*(subexp1);
end

function main()
    # simulation parameters
    scale = 2;                # set simulation size
    NX = 32*scale;            # domain size
    NY = NX;
    NSTEPS = 200*scale*scale; # number of simulation time steps
    NMSG   =  50*scale*scale; #show messages every NMSG time steps
    vis = false;              # show visualization; set to false for performance measurements
    NVIS = NMSG;              # show visualization every NVIS steps
    tau = 1;                  # relaxation time
    u_max =  0.04/scale;      # maximum velocity
    nu   =  (2*tau-1)/6;      # kinematic shear viscosity
    rho0 = 1;                 # rest density
    Re = NX*u_max/nu;         # Reynolds number; not used in the simulation itself

    # Lattice parameters; note zero direction is last
    w  = [1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36 4/9]; # weights
    cx = [1 0 -1  0 1 -1 -1  1 0]; # velocities, x components
    cy = [0 1  0 -1 1  1 -1 -1 0]; # velocities, y components

    x = (1:NX) .-0.5;
    y = (1:NY) .-0.5;
    # [X,Y] = meshgrid(x,y);
    X = [j for i in y, j in x]
    Y = [i for i in y, j in x]

    # Initialize populations
    (rho,u,v) = taylorgreen(0,X,Y,NX,NY,nu,rho0,u_max);

    f = zeros(size(rho, 1), size(rho, 2), 9);
    equilibrium!(f, rho, u, v);

    # initialize temporary variable
    fprop = zeros(size(f)); # set up propagation
    feq = zeros(size(f));

    # fprintf('Simulating Taylor-Green vortex decay\n');
    # fprintf('      domain size: #ux#u\n',NX,NY);
    # fprintf('               nu: #g\n',nu);
    # fprintf('              tau: #g\n',tau);
    # fprintf('            u_max: #g\n',u_max);
    # fprintf('             rho0: #g\n',rho0);
    # fprintf('        timesteps: #u\n',NSTEPS);
    # fprintf('       plot every: #u\n',NVIS);
    # fprintf('    message every: #u\n',NMSG);
    # fprintf('\n');

    E = rho.*(u.*u + v.*v);
    E = sum(E[:]);
    @show E
    # fprintf('#u,#g,#g,#g,#g\n',0,E,0,0,0);

    # tstart = tic;

    # main loop
    for t=1:NSTEPS
        # Periodic streaming of whole domain
        # see 'doc circshift' for info

        fprop[:,:,1] = circshift(f[:,:,1],[ 0  1]);
        fprop[:,:,2] = circshift(f[:,:,2],[ 1  0]);
        fprop[:,:,3] = circshift(f[:,:,3],[ 0 -1]);
        fprop[:,:,4] = circshift(f[:,:,4],[-1  0]);

        fprop[:,:,5] = circshift(f[:,:,5],[ 1  1]);
        fprop[:,:,6] = circshift(f[:,:,6],[ 1 -1]);
        fprop[:,:,7] = circshift(f[:,:,7],[-1 -1]);
        fprop[:,:,8] = circshift(f[:,:,8],[-1  1]);

        fprop[:,:,9] = f[:,:,9];

        # Compute macroscopic quantities
        # density
        rho = sum(fprop, dims=3);

        # momentum components
        jx = sum(fprop[:,:,[1, 5, 8]], dims=3) - sum(fprop[:,:,[3, 6, 7]], dims=3);
        jy = sum(fprop[:,:,[2, 5, 6]], dims=3) - sum(fprop[:,:,[4, 7, 8]], dims=3);

        # velocity components
        u = jx ./ rho;
        v = jy ./ rho;

        # Compute equilibrium distribution
        equilibrium!(feq, rho,u,v);

        # Collision step
        f = (1-1/tau)*fprop + (1/tau)*feq;

        if mod(t,NMSG) == 0
            # Calculate analytical solution
            (rhoa,uxa,uya) = taylorgreen(t,X,Y,NX,NY,nu,rho0,u_max);

            # Kinetic energy
            E = rho.*(u.*u + v.*v);
            E = sum(E[:]);
            @show E

            # Sum square errors
            rhoe2 = (rho-rhoa).*(rho-rhoa);
            sumrhoe2 = sum(rhoe2[:]);

            uxe2 = (u-uxa).*(u-uxa);
            sumuxe2 = sum(uxe2[:]);

            uye2 = (v-uya).*(v-uya);
            sumuye2 = sum(uye2[:]);

            rhoa2 = (rhoa .- rho0).*(rhoa .- rho0);
            sumrhoa2 = sum(rhoa2[:]);
            sumuxa2 = sum(uxa[:].*uxa[:]);
            sumuya2 = sum(uya[:].*uya[:]);

            # L2 norms
            L2rho = sqrt(sumrhoe2/sumrhoa2);
            L2ux = sqrt(sumuxe2/sumuxa2);
            L2uy = sqrt(sumuye2/sumuya2);

            # fprintf('#u,#g,#g,#g,#g\n',t,E,L2rho,L2ux,L2uy);
        end

        if vis && mod(t,NVIS) == 0
            figure(1)
            clf

            umag = sqrt(u.*u .+ v.*v) ./ u_max;
            imagesc(x ./ NX, y ./ NY, umag)
            # caxis([0 1])
            # colormap(hot(1024))
            hcb = colorbar;
            ylabel(hcb,'|u|/u_{max}')

            hstr = streamslice(X/NX,Y/NY,u,v);
            set(hstr,'Color','white');

            axis xy equal tight
            xlabel('x/l_x')
            ylabel('y/l_y')


            td = 1/(nu*(2*pi/NX)^2 + (2*pi/NY)^2);
            title(['flow field at t/t_d = ' num2str(t/td)])
            drawnow
        end
    end

    # Calculate performance information after the simulation is finished
    # runtime = toc(tstart);
    # nodes_updated = NSTEPS*NX*NY;
    # speed = nodes_updated/(1e6*runtime);

    # fprintf(' ----- performance information -----\n');
    # fprintf('        timesteps: #u\n',NSTEPS);
    # fprintf('          runtime: #.3f (s)\n',runtime);
    # fprintf('            speed: #.2f (Mlups)\n',speed);
end
