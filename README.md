# DropS: DROPlet Simulation

Author: Nikolaos Chamakos (nikoscham@gmail.com)
Institution: National Technical University of Athens, Greece

# Introduction

DropS is a Fortran-based simulation code designed to predict the equilibrium shape of a droplet on a structured substrate. It solves the Young-Laplace equation of capillary hydrostatics augmented with a disjoining pressure term, which models the solid-liquid interactions. This approach treats the liquid-solid and liquid-vapor interfaces in a unified framework, making it particularly effective for modeling cases where the number of contact lines is a priori unknown.

# Features

- Solvers:

  - Eikonal Equation Solver
  - Augmented Young-Laplace Equation Solver (Zero-Order Continuation)
  - Augmented Young-Laplace Equation Solver (Pseudo Arc-Length Continuation)
  - Decoupled Electrowetting Solver (Zero-Order Continuation)
  - Decoupled Electrowetting Solver (Pseudo Arc-Length Continuation)
  - Time Step Electrowetting Solver
  - Time Step Augmented Young-Laplace Equation Solver
  - Conventional Young-Laplace Equation Solver

- Documentation:
  - Chamakos, Kavousanakis, Papathanasiou, Soft Matter (2013)
  - Chamakos, Kavousanakis, Papathanasiou, Langmuir (2014)
  - Kavousanakis, Chamakos, Papathanasiou, J. Phys. Chem. C (2015)

# Installation

1. Prerequisites:

   - Fortran compiler (e.g., gfortran)
   - MUMPS library

2. Clone the repository:
   git clone <repository_url>
   cd DropS

3. Compile the code:
   make

# Usage

1. Prepare the input files:

   - GlobalParameters.txt
   - AYLParameters.txt
   - EIParameters.txt
   - CYLParameters.txt
   - ELParameters.txt

2. Run the simulation:
   ./DropS

3. View the results:
   - Results are saved in the Results directory.

# File Descriptions

Source Files:

- DropS.f90: Main program file.
- cvar.f90: Common variables module.
- auxfunct.f90: Auxiliary functions.
- initsol.f90: Initialization and solution reading routines.
- nodnum.f90: Nodal numbering routines.
- coord.f90: Coordinate transformation routines.
- defpar.f90: Parameter definition routines.
- tsfun.f90: Basis function evaluation routines.
- eikon.solv.f90: Eikonal equation solver.
- printres.f90: Result printing routines.
- ew.solv.f90: Electrowetting solver.
- ew.paramsolv.f90: Electrowetting solver with parametric continuation.
- ayl.paramsolv.f90: Augmented Young-Laplace solver with parametric continuation.
- ayl.solv.f90: Augmented Young-Laplace solver.
- cyl.solv.f90: Conventional Young-Laplace solver.

Input Files:

- GlobalParameters.txt: Global parameters for the simulation.
- AYLParameters.txt: Parameters for the Augmented Young-Laplace equation.
- EIParameters.txt: Parameters for the Eikonal equation.
- CYLParameters.txt: Parameters for the Conventional Young-Laplace equation.
- ELParameters.txt: Parameters for the Electrowetting simulation.

# License

DropS is distributed under the terms of the MIT license. See LICENSE for more information.

# Contact

For any questions or issues, please contact Nikolaos Chamakos at nikoscham@gmail.com.
