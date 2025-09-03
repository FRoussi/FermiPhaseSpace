# Fermionic Phase Space Representation Dynamics

The FermiPhaseSpace code computes solutions of the time-dependent Schrödinger equation on a Fermi-Hubbard system.

## Installation

I am not confident in this part. Please ask your favorite LLM.

And do not forget Pascal's wager : it costs you nothing to be polite to the machine, but might save you in the event of a robot uprising.

### Requirements

- Armadillo
- OPENBLAS
- ARPACK
- openMP

A MPI version of the fTWA method is avalaible.

We are using the [Armadillo library](https://arma.sourceforge.net/) for the linear algebra calculations. Its installation requires openBLAS and ARPACK, we have used and tested Armadillo with openBLAS as LAPACK.

### Programs

This code contains five programs, the two important ones are:

- "makeSitesSystem.cpp" that creates the entry files that define the Fermi-Hubbard systems.
- "fps.cpp", the main program that computes the dynamic of a Fermi-Hubbard model.

The others are:

- "fps_HF.cpp" that computes the dynamic of a Fermi-Hubbard model only with the Hartree-Fock method.
- "fps_fTWA_MPI.cpp", a MPI version of the Truncated Wigner Approximation dynamic computation.
- "makeWagnierFile.cpp" that turns a json Hartree-Fock system (the program entry) into a file usable for the Wanier tool suite.

### Makefile examples

The Makefile used to compile the four main executables are in the /FolderMakefile. They were used on a Cray environment. To compile the fps code, copy Makefile_fps into Makefile, modify the core directories, and run

    make

#### Options used

The options we used are not the consequences of a deep knowledge of computer sciences, but they worked well :

- Armadillo installation: '-DOPENBLAS_PROVIDES_LAPACK=true'
- FPS compilation : -larmadillo -fopenmp
- Environment variables : 'OPENBLAS_NUM_THREADS=1'

### Test

Once "makeSiteSytem.cpp"  and "fps.cpp" have been compiled, "tests.py" will generate a little Fermi-Hubbard system and try the four different methods on it.

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

- Cubic grid : -*grid* nx ny nz a-
- Bulk Centered Cubic : -*bcc* nx ny nz a-
- Face Centered Cubic : -*fcc* nx ny nz a-
- Tetragones : -*triangular* nx ny nz a-
- Graphene sheet build by Hexagones : -*honeyCombs* nc a-
- Graphene sheet build with the minimal cell of two sites : -*graphene* nx ny a-
- 2 graphene sheets with a phi angle : -*moire* nx ny phi d a-
- Kagome lattice : -*kagome* nx ny a-

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

We first make a system of 6 sites organised as a 3x2 grid :

    ./makeSitesSystem 6sites grid 3 2

We have now in the folder *systemsFolder* a two files called *6sites.json* and *6sites.xyz* that describes respectively the dynamical simulation and the sites positions. To launch the simulation :

    ./fps systemsFolder/6sites.json

This will compute a Hartree-Fock dynamic (default settings) and write a result file called *6sites_HF_res.json* in the same folder than *6sites.json*.
To modify the simulation settings, open *6sites.json*

## Fermions Phase-Space Representation methods

The FermiPhaseSpace code computes solutions of the time-dependent Schrödinger equation on a Fermi-Hubbard system.

- $` \frac{\partial}{\partial t} \ket{\psi} = -\frac{i}{\hbar} \hat{H} \ket{\psi}`$.

With $`  H `$ the Fermi-Hubbard Hamilltonian, written with number operators $` \hat{n}_{ij\sigma} = \hat{c}_{j\sigma}^{\dagger} \hat{c}_{i\sigma} `$ :

- $`  \hat{H}(t) = \sum_{i,j,\sigma} t_{ij} \hat{n}_{ij\sigma} + \sum_i u_i \hat{n}_{ii\sigma} \hat{n}_{ii\overline{\sigma}} + \sum_{i,j,\sigma} V_{ij}(t) \hat{n}_{ij\sigma}  `$.

We have implemented four different methods to solve the fermions dynamic :

### Exact Diagonalisation (ED)

The whole Hilbert space at constant fermion number is modeled. The exact dynamic is then written as $` \ket{\psi (t)} =  e^{-\frac{i}{\hbar}Ht} \ket{\psi_0} `$.

Using the Krylov subspace generated by $`H`$ and $` \ket{\psi} `$, we can approximate very accurately  $`e^{-\frac{i}{\hbar}Ht} \ket{\psi} `$ when $` t `$ is small. Hence we compute incrementally the ED dynamic with :

$` \ket{\psi (t + d t)} =  e^{-\frac{i}{\hbar}H d t} \ket{\psi (t)} `$.

Representing the whole Hilbert space is indeed costly. for half-hilled systems, the basis size scales exponentially with the number of sites : as $` O(4^n) `$. Depending on the computer's memory size and on the connectivity of the system, the method runs out of memory for systems of around 13-15 sites.

### Mean-field, Hartree-Fock (HF)

The usual simplification is to consider the one-particle density matrix $` \rho_{ij\sigma} =  \braket{ \hat{n}_{ij\sigma}} `$ and compute the Coulomb interaction as the interaction between the particle and the electric potential it creates. The Hamiltonian in its one-particle form is

$` H_{ij} =  t_{ij} + u_i \rho_{ij\overline{\sigma}} `$

The time derivative is computed from the Heisenberg picture

$` \frac{d}{dt} \rho =  \frac{i}{\hbar} [H, \rho] `$

and the density matrix trajectory is with a Runge-Kutta algorithm.

### Fermionic Truncated Wigner Approximation (fTWA)

The phase-space variables are mapped to the anti-symmetrised number operators:

$` \rho_{ij\sigma} <->  \braket{\hat{E}_{i\sigma}^{j\sigma}} = \frac{1}{2} \braket{\hat{c}^{\dagger}_{j\sigma}\hat{c}_{i\sigma} - \hat{c}_{i\sigma}\hat{c}^{\dagger}_{j\sigma}} `$

The wave-function is mapped to a density probability of one-particle density matrix, and the time-dependent Schrödinger equation becomes a density partial differential equation without even order terms. Hence the truncation at the third order gives us a linear differential equation very close to the mean-field approximation. The density is modeled by with cloud of points, each one of them following independently its own trajectory. Averages of this cloud of points gives us the usual observable like energy, site occupation and site correlation.

### Gaussian Phase-Space Representation (GPSR)

This time, the phase-space variable are mapped to normally-ordered number operators:

$` \rho_{ij\sigma} <->  \braket{\hat{c}^{\dagger}_{j\sigma}\hat{c}_{i\sigma}} `$

The terms above second order of the density partial differential equation disappear so we find a Fokker-Planck equation. This can be modeled by trajectories following stochastic differential equations.

$` d \rho_{ij\sigma} =  A_{ij\sigma}(\rho) dt + B_{ij\sigma}(\rho) dW_t `$

## Acknowledgements

This code was developped the supervision Olle Eriksson, Uppsala University, with the support from the Swedish Research Council, the Knut and Alice Wallenberg Foundation, the European Research Council via Synergy Grant No. 854843- FASTCORR, Swedish energy agency and eSSENCE.

---- ***Capitaine Docteur Frousse***
