# HDFIT HPC Applications Suite for OpenBLAS

This directory contains scripts to clone, patch and build the HPC applications supported by HDFIT, as well as prepare the input configurations available for each application. The selected applications use OpenBLAS and are suitable for GEMM fault injection experiments.

## Available Applications

As of July 5th 2022, the following applications are supplied with HDFIT:

* __QMCPack__: a many-body ab initio Quantum Monte Carlo code for computing the electronic structure of atoms, molecules, 2D nanomaterials and solids. It can be compiled with the __qmcpack__ make target.
* __CP2K__: a quantum chemistry and solid-state physics software package that can perform atomistic simulations of solid state, liquid and molecular systems. It can be compiled with the __cp2k__ make target.
* __NWChem__: a quantum chemistry package for biomolecules, nano structures and solid-state, supporting classical and quantum simulations, plane-wave or Gaussian basis functions. It can be compiled with the __nwchem__ make target.
* __Quantum Espresso__: a suite of codes for electronic-structure calculations and materials modeling at the nanoscale, based on density-functional theory. It can be compiled with the __qe__ make target.
* __Remhos__: a mini-app that solves the pure advection equations that are used as part of the Eulerian phase in Arbitrary Lagrangian-Eulerian (ALE) simulations. It can be compiled with the __remhos__ make target.
* __High-Performance Linpack__: a benchmark which solves a random dense linear system in double precision arithmetic, on distributed-memory computers. It can be compiled with the __hpl__ make target.

## External Dependencies

Compiling the HPC applications requires a basic set of generic dependencies. These are __make__, __cmake__, __autoconf__, __pkgconf__ and a functional __gcc__ and __gfortran__ toolchain. Furthermore, a valid __MPI__ installation is required (OpenMPI is recommended). The specific additional dependencies for each HPC application are compiled automatically as part of the build process, and no user action is required. A functional __Python 3__ installation is recommended.

Before attempting to compile the HPC applications, make sure to configure the correct path to the OpenBLAS HDFIT library (__OPENBLAS\_ROOT__ variable) 
within the __config.mk__ file, as indicated in the main [README](../README.md) document.

## Compilation Steps

In order to build all HPC applications, the following command can be used:

```
make all
```

Building specific HPC applications can be achieved by using the corresponding make targets. This repository does not contain the complete source code for each application, but rather only the patches necessary to enable HDFIT use: as such, the build process first clones each application's source code from the respective repositories, then applies the HDFIT patches, and finally compiles them. Any dependencies are compiled as part of this process.

Please note that it may be necessary to use the __LD\_PRELOAD__ environment variable (set to the path to the HDFIT OpenBLAS .so library file) for HPC applications to work correctly. This is to prevent using a different version of OpenBLAS already present on the system, as well as to address the specific linking mechanisms in each HPC application. In case errors are encountered during the build process, please refer to the README documents in each application's directory for more detailed compilation instructions, and to the respective documentation websites.

The supplied HPC applications have been compiled and tested on Ubuntu (18.04, 20.04 and 22.04), as well as on SLES (12SP5). Compatibility with Windows WSL has been verified as well.

# Application Configurations

In the following is an overview of the available configurations supplied with HDFIT for each application, as well as the external data dependencies
they have, if any. It should be noted that each supplied application configuration comes with an associated HDFIT testing script, using the naming
convention __APPLICATION-test-CONFIGURATION.env__, which can be used as input to the __HDFIT\_runner.sh__ script to perform fault injection experiments.
Please refer to the README documents of each separate application for additional information.

## QMCPack

QMCPack's dependencies are LibXML, FFTW, HDF5 and Boost (headers included), which are installed automatically by HDFIT. The following two
configurations are supplied:

* __in.NiO8__: Nickel Oxide with 32 atoms (__NiO-fcc-S8-dmc.xml__ file, requires __NiO-fcc-supertwist111-supershift000-S8.h5__);
* __in.NiO16__: Nickel Oxide with 64 atoms (__NiO-fcc-S8-dmc.xml__ file, requires __NiO-fcc-supertwist111-supershift000-S16.h5__).

Please note that the extra files required for the Nickel Oxide simulations, which are to be placed in the same directory as the respective .xml
input files, can be found at [this](https://anl.app.box.com/s/pveyyzrc2wuvg5tmxjzzwxeo561vh3r0) link. Here is an example for the execution of a
Nickel Oxide simulation using 64 atoms:

```
cd in.NiO16
mkdir out && cd out
OMP_NUM_THREADS=1 mpirun -np 1 ../../build/bin/qmcpack ../NiO-fcc-S16-dmc.xml
```

## CP2K

CP2K requires the Libint, Libxc, Libxsmm, HDF5 and FFTW libraries, among others. The toolchain script included with CP2K will install all dependencies
automatically, while HDFIT will handle other remaining dependencies (such as the ZLib library). Our distribution of the application comes with the
two following sample configurations:

* __in.H2O__: Water molecule, energy and forces calculation (__H2O.inp__ file);
* __in.C2H4__: Ethene molecule, energy and forces calculation (__C2H4.inp__ file).

No molecular dynamics inputs were included, since the underlying algorithms make very little use of GEMM operations. Here is an example for the
execution of an Ethene simulation:

```
cd in.C2H4
mkdir out && cd out
OMP_NUM_THREADS=1 mpirun -np 1 ../../exe/local/cp2k.popt -i ../C2H4.inp
```

## NWChem

Aside from MPI, NWChem does not have any specific additional dependencies. We supply the following two sample configurations for the application:

* __in.Na16__: 16-atom Sodium Car-Parrinello simulation (__Na16.nw__ file);
* __in.3carbo__: 3-Carboxybenzisoxazole quantum molecular dynamics simulation (__3carbo\_dft.nw__ file).

Here is an example for the execution of a 3-Carbo simulation:

```
cd in.3carbo
mkdir out && cd out
OMP_NUM_THREADS=1 mpirun -np 1 ../../bin/LINUX64/nwchem ../3carbo_dft.nw
```

## Quantum Espresso

QE requires the FFTW library - aside from this, the application installs its other dependencies automatically, hence no further configuration is needed.
This distribution comes with two different sample configurations:

* __in.SrVO3__: Electronic structure calculation for SrVO3 Perovskite (__srvo3.scf.in__ file);
* __in.AlAs__: Electronic structure calculation for Aluminium Arsenide (__alas.scf.efield2.in__ file).

Please note that both simulations require extra pseudopotentials files that need to be downloaded externally. This can be performed automatically
by running the __run\_example__ script within the __PW/examples/example14__ and __PW/examples/example10__ directories respectively. After performing
this, application runs can be performed - here is an example for the execution of an AlAs simulation:

```
cd in.AlAs
mkdir out && cd out
OMP_NUM_THREADS=1 mpirun -np 1 ../../bin/pw.x -i ../alas.scf.efield2.in
```

## Remhos

Remhos requires the MFEM, Metis and Hypre libraries, which are installed automatically by HDFIT. This distribution comes with two different configurations:

* __in.remap__: Performs a simulation in Remap mode;
* __in.transport__: Performs a simulation in Transport mode.

Remohos does not require any external input files, and runs can be configured through command-line arguments. The settings chosen for the two
configurations can be found in the respective HDFIT testing scripts, __REMHOS-test-remap.env__ and __REMHOS-test-transport.env__. Here is an example
for the execution of a remap simulation:

```
cd in.remap
mkdir out && cd out
OMP_NUM_THREADS=1 mpirun -np 1 ../../remhos -m ../../data/cube01_hex.mesh -p 10 -rs 1 -rp 1 -o 2 -dt 0.02 -tf 0.8 -ho 1 -lo 2 -fct 2 -visit -vs 5
```

## High-Performance Linpack

Save for MPI, HPL does not have any additional dependencies that need to be handled by HDFIT. This distribution comes with two different configurations:

* __in.small__: A small configuration (__HPL.dat__ file);
* __in.large__: A larger configuration (__HPL.dat__ file).

Here is an example for the execution of an HPL run under the small configuration:

```
cd in.small
OMP_NUM_THREADS=1 mpirun -np 1 ../testing/xhpl HPL.dat
```

