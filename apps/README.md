# HDFIT HPC Applications Suite

This directory contains scripts to clone, patch and build the HPC applications supported by HDFIT, as well as prepare the input configurations available for each application.

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

Obviously, the HPC applications need to be linked to the custom HDFIT OpenBLAS library. To do so, replace the __PATH\_TO\_CUSTOM\_OPENBLAS\_LIB__ string in __config.mk__ with the path to the OpenBLAS root directory on your system, before proceeding with the build process.

## Compilation Steps

In order to build all HPC applications, the following command can be used:

```
make all
```

Building specific HPC applications can be achieved by using the corresponding make targets. This repository does not contain the complete source code for each application, but rather only the patches necessary to enable HDFIT use: as such, the build process first clones each application's source code from the respective repositories, then applies the HDFIT patches, and finally compiles them. Any dependencies are compiled as part of this process.

Please note that it may be necessary to use the __LD\_PRELOAD__ environment variable (set to the path to the HDFIT OpenBLAS .so library file) for HPC applications to work correctly. This is to prevent using a different version of OpenBLAS already present on the system, as well as to address the specific linking mechanisms in each HPC application. In case errors are encountered during the build process, please refer to the README documents in each application's directory for more detailed compilation instructions, and to the respective documentation websites.

The supplied HPC applications have been compiled and tested on Ubuntu (18.04, 20.04 and 22.04), as well as on SLES (12SP5). Compatibility with Windows WSL has been verified as well.

## Additional Data Requirements

In order to be used with the supplied input configurations, some of the HPC applications need additional data files from external sources. This is required for the following applications:

* __QMCPack__: for each input configuration, an HDF5 data file must be downloaded and placed alongside the corresponding .xml input file (i.e., in the __in.NiO8__ and __in.NiO16__ directories).
* __Quantum Espresso__: a series of pseudopotentials files must be downloaded into the __pseudo__ sub-directory.

Detailed instructions on how to retrieve the data files mentioned above can be found in each HPC application's README file, under the __Running Experiments__ section.
