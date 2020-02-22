# LatticeBoltzmann.jl

| **Build Status**                                                                                            | **Documentation**                 |
|:------------------------------------------------------------------------------------------------------------|:----------------------------------|
| ![Lifecycle][lifecycle-img] [![CI][ci-status-img]][ci-status-url] [![codecov.io][codecov-img]][codecov-url] | [![][docs-dev-img]][docs-dev-url] |

This package was developed for my thesis *A Theoretical & Practical Framework for Lattice Boltzmann Models*. Currently the documentation is a work in progress.

The [**notebooks**](https://github.com/MarkRedeman/LatticeBoltzmann.jl/tree/master/examples/notebooks) contain reproducible examples of results that I've used in my thesis.
I've focused on four problems:
- [**Shear wave decay**](https://github.com/MarkRedeman/LatticeBoltzmann.jl/blob/master/examples/notebooks/shear_wave.ipynb): shows we solve initial value problems and can obtain steady state solutions with quadratic convergence.
- [**Taylor Green Vortex**](https://github.com/MarkRedeman/LatticeBoltzmann.jl/blob/master/examples/notebooks/taylor_green_vortex.ipynb): Shows we can obtain quadratic convergence for the velocity and stress components. It also shows the effect the initialization strategy of a Lattice Boltzmann Method can have on its accuracy (in particular it shows that having an inconsistent pressure can lead to larger errors).
- [**Couette flow**](https://github.com/MarkRedeman/LatticeBoltzmann.jl/blob/master/examples/notebooks/couette.ipynb): Shows that the D2Q9 quadrature can find an exact solution of the Couette flow. The other quadratures are only second or first order accurate due to an inexact boundary condition.
- [**Poiseuille**](https://github.com/MarkRedeman/LatticeBoltzmann.jl/blob/master/examples/notebooks/poiseuille.ipynb): shows how we can find a solution to the 2D poiseuille flow and how its accuracy depends on the relaxation time. 


## Quadratures
- D2Q4 
- D2Q5 
- D2Q9 
- D2Q13
- D2Q17
- D2Q21
- D2Q37

## Collision models
- SRT (+ Force)
- TRT (+ Force)
- Regularized (projected) MRT (slow) (+ Force)

## Initial conditions
- Velocity (equilibrium w. constant density)
- Velocity + stress (offequilibrium)
- Velocity + pressure (analytical equilibrium)
- Velocity + pressure + stress (analytical equilibrium + offequilibrium)
- Iterative (Mei et al)

## Boundary conditions
- Periodic (default)
- Bounce back
- Moving wall

## Test problems
- Shear wave (decay / steady state)
- Taylor Green Vortex (decay / steady state)
- Couette
- Poiseuille

## TODO

### Project
- Add documentation & doctests
- Improve code coverage
- Try [mutation testing](https://github.com/MikeInnes/Vimes.jl)

### Features
- Add abstraction for simulation parameters (Δx, Δy, Re)
- Double Distribution Functions for temperature
- Double Distribution Functions for multicomponent and multiphase flows
- Improve plotting and process managers (the ugly code)
- Allow using different equilibrium functions
- Implement 1D and 3D features

#### Boundary conditions
Boundary conditions should be rewritten so that they can be applied to a node at
a location not necessarily at the boundary, so that we can have second order
bounce back boundary conditions.

- Pressure difference
- Temperature
- Multispeed boundary conditions

### Performance / infrastructure
- Improve benchmarks
- Refactor to use [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl/)
- Refactor collision models to use kernel methods (will make it easier to allow distributed computing)
- Add instructions for Docker and Singularity

 
<!-- References and urls -->
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://markredeman.github.io/LatticeBoltzmann.jl/

[ci-status-url]: https://github.com/MarkRedeman/LatticeBoltzmann.jl/actions?query=workflow%3ACI
[ci-status-img]: https://github.com/MarkRedeman/LatticeBoltzmann.jl/workflows/CI/badge.svg

[codecov-img]: http://codecov.io/github/MarkRedeman/LatticeBoltzmann.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/MarkRedeman/LatticeBoltzmann.jl?branch=master

[lifecycle-img]: https://img.shields.io/badge/lifecycle-experimental-orange.svg

<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
