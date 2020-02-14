# LatticeBoltzmann.jl

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![Build Status](https://travis-ci.com/MarkRedeman/LatticeBoltzmann.jl.svg?branch=master)](https://travis-ci.com/MarkRedeman/LatticeBoltzmann.jl)
[![codecov.io](http://codecov.io/github/MarkRedeman/LatticeBoltzmann.jl/coverage.svg?branch=master)](http://codecov.io/github/MarkRedeman/LatticeBoltzmann.jl?branch=master)

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
- Double Distribution Functions for temperature
- Double Distribution Functions for multicomponent and multiphase flows
- Improve plotting and process managers (the ugly code)
- Allow using different equilibrium functions
- Implement 1D and 3D features

### Performance / infrastructure
- Improve benchmarks
- Refactor to use [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl/)
- Refactor collision models to use kernel methods (will make it easier to allow distributed computing)
- Add instructions for Docker and Singularity
