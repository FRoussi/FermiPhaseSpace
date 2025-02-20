# THAT IS MY COOL PROJECT

The FermiPhaseSpace code computes solutions of the time-dependent Schrödinger equation on a Fermi-Hubbard system.

## Installation

I am not confident in this part. Please ask your favorite LLM.

And do not forget Pascal's bet : it costs you nothing to be polite to the machine, but might save you in the event of a robot uprising.

### Requirements

- Armadillo
- OPENBLAS
- ARPACK
- openMP

A MPI version of the fTWA method is avalaible.

### Compilation of programs

This code contains five programs whose main functions are in "fps.cpp", "fps_HF.cpp", "fps_fTWA_MPI.cpp", "makeSitesSystem,cpp" and "makeWannierFile.cpp".

### Makefile examples

The Makefile used to compile the four main executables are in the /FolderMakefile. They were used on a Cray environment.

## Usage

The workflow is supposed to be :

1. Create a Fermi-Hubbard system, described by a .json and a .xyz files.
2. Modify the .json file at your will as it contains all the simulation parameters (method, energy constants, step times etc ...).
3. Solve the dynamic.

### The FPS Fermi-Hubbard system

Fermi-Hubbard simulations are entirely decribed by a .json file, this file is the sole entry for the "fps" program and contains everything necessary. The json entries are :

- "sim_name"
- "method"
- "system"
- "parameters"
- "results"

#### Simulation name

The "sim_name" json entry is the name of the simulation. It is used to write the result file, called *sim_name*_*method*_res.json.

#### Simulation's methods

The "method" can be one of the four dunamical method implemented :

- "HF" : Hartree-Fock
- "fTWA" : Fermionic Truncated Wigner Approximation
- "GPSR" : Gaussian Phase-Space Representation
- "ED" : Exact Diagonalisation

#### Parameters

They are the parameters of the dynamical simulation.

- "dt" : time step between two snapshot of the simulation.
- "t_max" : time extension of the simulation.
- "stepsPerTimePoints" : number of simulation steps between two snapshots of the simulation, so dt/stepsPerTimePoints is the real simulation time step.
- "int_order" : for HF and fTWA methods, order of the Runge-Kutta integrator. For GPSR, order of the Runge-Kutta integrator of the linear part of its time steps.
- "n_traj" : for GPSR and fTWA, number of trajectory used to model the probability density.

#### Results

To increase efficiency, all possible observable are not computed. Hence we choose beforehand which variable we want to extract of the simulation. The "results" json entry is a list of boolean  telling if an observable should be computed and written in the result file.

- "diagDMatrix"
- "diag2Matrix"
- "onSiteDoublons"
- "wholeDMatrix"
- "energies"
- "mfStates"
- "tbStates"

#### System

The "system" json entry describes the Fermi-Hubbard model : its geometry, its energy functions and the simulation's initial condition.

- "lattice"
- "coulombInteraction"
- "tunnelingInteraction"
- "externalPotential"
- "IC"

### makeSitesSystem

The first program is "./makeSitesSystem". It takes as entry the name of the new system, the geometry of the sites lattice, and several parameters for this geometry. If the geometry is recognised, two files are created in a folder called "systemsFolder". The .json file of the simulation with generic parameters, and the .xyz file with the site positions. Both files bear the name given as parameter. The use should be :

    ./makeSitesSystem -nameOfTheSystem- -geometryOfTheSites- -param1- -param2- ...

For example, "./makeSitesSystem 6sites grid 3 2" writes two files called "6sites.json" and "6sites.xyz" describing a 3x2 cubic system.

The possible lattice geometry are

- grid nx ny nz a               cubic grid
- bcc nx ny nz a                Bulk Centered Cubic
- fcc nx ny nz a                Face Centered Cubic
- triangular nx ny nz a         Tetragones
- honeyCombs nc a               Graphene sheet build by Hexagones
- graphene nx ny a              Graphene sheet build with the minimal cell of two sites
- moire nx ny phi d a           2 graphene sheets with a phi angle
- kagome nx ny a                Kagome lattice

Missing parameters are replaced by "1".

### fps

    ./fps -jsonFile-

### fps_HF

    ./fps_HF -jsonFile-

fps_HF will ignore the json entry "method" and run a Hartree-Fock simulation.

### fps_fTWA_MPI

    mpirun -n -nProc- ./fps_fTWA -jsonFile-

fps_HF will ignore the json entry "method" and run a fTWA simulation on -nProc- processus using MPI.

### makeWannierFile

Is used to generate files for the Wannier tools suite.

### Example

## Fermions Phase-Space Representation methods

The FermiPhaseSpace code computes solutions of the time-dependent Schrödinger equation on a Fermi-Hubbard system.

- $` \frac{\partial}{\partial t} \ket{\psi} = -\frac{i}{\hbar} \hat{H} \ket{\psi}`$.

With $`  H `$ the Fermi-Hubbard Hamilltonian, written with number operators $` \hat{n}_{ij\sigma} = \hat{c}_{j\sigma}^{\dagger} \hat{c}_{i\sigma} `$ :

- $`  \hat{H}(t) = \sum_{i,j,\sigma} t_{ij} \hat{n}_{ij\sigma} + \sum_i u_i \hat{n}_{ii\sigma} \hat{n}_{ii\overline{\sigma}} + \sum_{i,j,\sigma} V_{ij}(t) \hat{n}_{ij\sigma}  `$.

We have implemented four different methods to solve the fermions dynamic :

### Exact Diagonalisation (ED)

The whole Hilbert space at constant fermion number is modeled. The exact dynamic is then written as $` \ket{\psi (t)} =  e^{-\frac{i}{\hbar}Ht} \ket{\psi_0} `$.

Using the Krylov subspace generated by $`H`$ and $` \ket{\psi} `$, we can approximate very accurately  $`e^{-\frac{i}{\hbar}Ht} \ket{\psi} `$ when $` t `$ is small. Hence we compute incrementally the ED dynamic with :

$` \ket{\psi (t + \partial t)} =  e^{-\frac{i}{\hbar}H \partial t} \ket{\psi (t)} `$.

Representing the whole Hilbert space is indeed costly. for half-hilled systems, the basis size scales exponentially with the number of sites : as $` O(4^n) `$. Depending on the computer's memory size and on the connectivity of the system, the method runs out of memory for systems of arounf 13-15 sites.

### Mean-field, Hartree-Fock (HF)

The usual simplification is to consider the one-particle density matrix $` \rho_{ij\sigma} =  < \hat{n}_{ij\sigma}> `$ and compute the Coulomb interaction as the interaction between the particle and the electric potential it creates.

### Fermionic Truncated Wigner Approximation (fTWA)

### Gaussian Phase-Space Representation (GPSR)

## Acknowledgements

Thanks ME

---- Capitaine Docteur Frousse
