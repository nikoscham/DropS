## DropS

DropS is a Fortran-based simulation code designed to predict the equilibrium shape of a droplet on a structured substrate. It solves the Young-Laplace equation of capillary hydrostatics augmented with a disjoining pressure term, which models the solid-liquid interactions. This approach treats the liquid-solid and liquid-vapor interfaces in a unified framework, making it particularly effective for modeling cases where the number of contact lines is a priori unknown.

## Documentation

DropS has been used in several scientific papers:
  - <a href="https://pubs.rsc.org/en/content/articlelanding/2013/sm/c3sm51377g" target="_blank">Chamakos, Kavousanakis, Papathanasiou, Soft Matter (2013)</a>
  - <a href="https://pubs.acs.org/doi/full/10.1021/la500408j" target="_blank">Chamakos, Kavousanakis, Papathanasiou, Langmuir (2014)</a>
  - <a href="https://pubs.acs.org/doi/abs/10.1021/acs.jpcc.5b00718" target="_blank">Kavousanakis, Chamakos, Papathanasiou, J. Phys. Chem. C (2015)</a>

## How to Use

1. Prepare the input files:

   - GlobalParameters.txt
   - AYLParameters.txt
   - EIParameters.txt
   - CYLParameters.txt
   - ELParameters.txt

2. Compile the code:
   ```sh
   make
   ```

3. Run the simulation:
   ```sh
   ./DropS
   ```

4. The results are saved in the Results directory

### Prerequisites

   - Fortran compiler: Install a Fortran compiler, such as <a href="https://gcc.gnu.org/fortran/" target="_blank">GFortran</a>
   - MUMPS Library: Ensure the <a href="https://mumps-solver.org/index.php" target="_blank">MUMPS</a> library  is installed in the directory <code>/usr/local/MUMPS</code>
   - Linux Environment: The application requires a Linux operating system for compatibility

## File Descriptions

Source Files:

- DropS.f90: Main program file
- cvar.f90: Common variables module
- auxfunct.f90: Auxiliary functions
- initsol.f90: Initialization and solution reading routines
- nodnum.f90: Nodal numbering routines
- coord.f90: Coordinate transformation routines
- defpar.f90: Parameter definition routines
- tsfun.f90: Basis function evaluation routines
- eikon.solv.f90: Eikonal equation solver
- printres.f90: Result printing routines
- ew.solv.f90: Electrowetting solver
- ew.paramsolv.f90: Electrowetting solver with parametric continuation
- ayl.paramsolv.f90: Augmented Young-Laplace solver with parametric continuation
- ayl.solv.f90: Augmented Young-Laplace solver
- cyl.solv.f90: Conventional Young-Laplace solver

Input Files:

- GlobalParameters.txt: Global parameters for the simulation
- AYLParameters.txt: Parameters for the Augmented Young-Laplace equation
- EIParameters.txt: Parameters for the Eikonal equation
- CYLParameters.txt: Parameters for the Conventional Young-Laplace equation
- ELParameters.txt: Parameters for the Electrowetting simulation

## License

DropS is distributed under the terms of <a href="./LICENSE" target="_blank">MIT License</a>. &#169; 2024 Nikolaos Chamakos.
