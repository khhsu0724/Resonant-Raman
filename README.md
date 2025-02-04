# Resonant Raman spectroscopy
This project provides a codebase designed to calculate **Resonant Raman Scattering** with Hubbard like Models
The repository contains 3 codes:
- ResRaman-hny.f90: An MPI code that can calculate around 18~20 sites Hubbard model on Square/Triangular/Honeycomb/Kagome lattices. The code is also capable of calculating Non-Resonant Raman spectroscopy and optical conductivity.
- Heisenberg.ipynb: A simple python code that calculates effective Resonant Raman scattering on the Heisenberg model using effective Raman scattering operator (see below).
- Hopping.ipynb: Calculates effective Raman scattering operator in the large U limit and off-resonance using [Shastry-Shraiman's technique](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.65.1068).


## Package Dependencies
- **fortran**
- **MPI**
- **BLAS/LAPACK**
- **ARPACK**

## Compile and running
The fortran code is only tested on NERSC
```
$ ftn -L$ARPACKROOT/lib -larpack -lparpack -llapack -O3 -fallow-argument-mismatch -fbounds-check -o ResRaman ResRaman-hny.f90
```
## Publication
arxiv post coming
