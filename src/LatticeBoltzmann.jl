module LatticeBoltzmann

import Base: range

using LinearAlgebra
using Plots

δ(α, β) = α == β ? 1 : 0

include("quadratures.jl")

# Equilibria
include("hermite_polynomials.jl")

# Velocity Distribution Function related methods
include("velocity_distribution_function/moments.jl")
include("velocity_distribution_function/maxwell_boltzmann.jl")
include("velocity_distribution_function/hermite.jl")
include("velocity_distribution_function/quadratures.jl")

include("boundary_conditions.jl")
include("problems/problems.jl")
include("initial_conditions.jl")
include("stream.jl")
include("collision_models.jl")
include("processing_methods.jl")
include("lattice_boltzmann_model.jl")

include("interop/plots.jl")

export CollisionModel, SRT, TRT, MRT
export Quadrature, D2Q4, D2Q5, D2Q9, D2Q13, D2Q17, D2Q21, D2Q37, opposite

export analayze_convergence
# Problems
export DecayingShearFlow
export LidDrivenCavityFlow
export TaylorGreenVortex
export CouetteFlow
export GenericFluidFlowProblem
export LinearizedThermalDiffusion, LinearizedTransverseShearWave
export PoiseuilleFlow

export process!,
    apply_boundary_conditions!,
    density,
    velocity,
    pressure,
    temperature,
    decay,
    force,
    FluidFlowProblem,
    viscosity,
    delta_t
export Lattice,
    initialize,
    AnalyticalEquilibriumAndOffEquilibrium,
    AnalyticalEquilibrium,
    AnalyticalVelocity,
    IterativeInitialization,
    # Thermodynamics
    density,
    momentum,
    pressure,
    total_energy,
    kinetic_energy,
    internal_energy,
    temperature,
    dimension,
    # Equilibria
    equilibrium,
    equilibrium!,
    hermite_equilibrium,
    hermite_first_nonequilibrium,
    hermite,
    # Siumlation
    stream!,
    collide!,
    simulate

end
