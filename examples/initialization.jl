using lbm, Plots, DataFrames

q = D2Q9()
τ = 1.0
problem = lbm.TGV(q, τ, 1//2, 2 * 8, 2 * 16)

τ = 0.8
problem = lbm.TGV(q, τ, 1//2, 2 * 16, 2 * 16)

# There still seems to be an scaling issue
# problem2 = TaylorGreenVortex(1.0 / 6.0, 1//2, static = false, a = 2);

iss = [

    # With this initialization strategy it is assumed that the initial density field
    # (or equivalently the pressure field p) is not available
    lbm.ConstantDensity(),

    # The initial pressure field p is known, which, using the equation fo state,
    # is used to set the initial density
    lbm.AnalyticalEquilibrium(),

    # Both an initial pressure field and gradient of the velocity is known.
    # The gradient of the velocity is used to initialize the offequilibrium components
    lbm.AnalyticalEquilibriumAndOffEquilibrium(),

    # We initialize f using an iterative procedure where only the density is conserved
    # it was shown that this procedure gives consistent itnitial conditions for both
    # the equilibrium and off equilibrium components
    lbm.IterativeInitializationMeiEtAl(τ, 1E-9),
]

result = lbm.simulate(problem, q, t_end = 100, initialization_strategy = iss[1])
