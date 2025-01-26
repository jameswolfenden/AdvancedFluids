# Fluid Dynamics Solver

A parallel implementation of a 3D fluid dynamics solver in C++. Outputs simulation data in HDF5 format for visualisation.

## Key Classes

1. **RiemannSolver**
   - Uses iterative Riemann solver, with Newton-Raphson method
   - Computes fluxes at cell interfaces

2. **DomainEulerSolver**
   - Works to update domains using the Riemann solver
   - Parallelises the domain update process

### I/O System

1. **HDF5Writer**
   - Writes simulation data to HDF5 files

2. **XDMFWriter**
   - Generates visualisation metadata
   - Links HDF5 data files

## Acknowledgements

Implementation based on the numerical methods described in:
> Toro, E. F. (2009). Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction.