using LatticeBoltzmann, Plots, DataFrames

function main()
    q = D2Q9()
    τ = 1.0
    problem = LatticeBoltzmann.TGV(q, τ, 1//2, 2 * 8, 2 * 16)

    τ = 0.8
    problem = LatticeBoltzmann.TGV(q, τ, 1//2, 2 * 16, 2 * 16)

    # There still seems to be an scaling issue
    # problem2 = TaylorGreenVortex(1.0 / 6.0, 1//2, static = false, a = 2);

    iss = [

        # With this initialization strategy it is assumed that the initial density field
        # (or equivalently the pressure field p) is not available
        LatticeBoltzmann.ConstantDensity(),

        # The initial pressure field p is known, which, using the equation fo state,
        # is used to set the initial density
        LatticeBoltzmann.AnalyticalEquilibrium(),

        # Both an initial pressure field and gradient of the velocity is known.
        # The gradient of the velocity is used to initialize the offequilibrium components
        LatticeBoltzmann.AnalyticalEquilibriumAndOffEquilibrium(),

        # We initialize f using an iterative procedure where only the density is conserved
        # it was shown that this procedure gives consistent itnitial conditions for both
        # the equilibrium and off equilibrium components
        LatticeBoltzmann.IterativeInitializationMeiEtAl(τ, 1E-10),
    ]

    result = LatticeBoltzmann.simulate(problem, q, t_end = 100, initialization_strategy = iss[1])
end

function initialization_problem()
    # Book:
    (Nx, Ny) = (96, 72)
    Δt = 1.0
    Δx = 2pi / Nx
    τ = 0.8 Δt

    # ν = (τ - 0.5) / cs = 0.3 Δt / cs = 0.1 Δt
    ν = 0.1 Δx^2 / Δt

    û₀ = 0.03 Δx / Δt
    ρ̄ = 1.0
    p₀ = 0.0
    t_d = 840 Δt


    problem = LatticeBoltzmann.TGV(q, τ, 2, Nx, Ny, û₀)
end

function main_book()
    iss = [
        # With this initialization strategy it is assumed that the initial density field
        # (or equivalently the pressure field p) is not available
        LatticeBoltzmann.ConstantDensity(),

        # The initial pressure field p is known, which, using the equation fo state,
        # is used to set the initial density
        LatticeBoltzmann.AnalyticalEquilibrium(),

        # Both an initial pressure field and gradient of the velocity is known.
        # The gradient of the velocity is used to initialize the offequilibrium components
        LatticeBoltzmann.AnalyticalEquilibriumAndOffEquilibrium(),

        # We initialize f using an iterative procedure where only the density is conserved
        # it was shown that this procedure gives consistent itnitial conditions for both
        # the equilibrium and off equilibrium components
        LatticeBoltzmann.IterativeInitializationMeiEtAl(τ, 1E-10),
    ]

    # Book:
    (Nx, Ny) = (96, 72)
    Δt = 1.0
    Δx = 2pi / Nx
    τ = 0.8 Δt

    # ν = (τ - 0.5) / cs = 0.3 Δt / cs = 0.1 Δt
    ν = 0.1 Δx^2 / Δt

    û₀ = 0.03 Δx / Δt
    ρ̄ = 1.0
    p₀ = 0.0
    t_d = 840 Δt


    problem = LatticeBoltzmann.TGV(q, τ, 2, Nx, Ny, û₀)

    Δt = delta_t(problem)
    t_end = 840 * Δt

    result = LatticeBoltzmann.simulate(problem, q, t_end = t_d, initialization_strategy = iss[1])
end
