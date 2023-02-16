# HDFIT + LLTFI HPC Applications Suite

This directory contains scripts to clone and build a set of HPC applications supported by HDFIT, as well as prepare the input configurations 
available for each of them. This sub-directory focuses on LLVM-based compiler instrumentation for FPU fault injection: this is achieved through integration between the HDFIT and [LLTFI](https://github.com/DependableSystemsLab/LLTFI) toolchains.

## Available Applications

As of November 11th 2022, the following applications are supplied with HDFIT:

* __QMCPack__: a many-body ab initio Quantum Monte Carlo code for computing the electronic structure of atoms, molecules, 2D nanomaterials and solids. It can be compiled with the __qmcpack__ make target.
* __GROMACS__: a molecular dynamics package mainly designed for simulations of proteins, lipids, and nucleic acids. It can be compiled with the __gromacs__ make target.
* __SeisSol__: an open-source software for the simulation of seismic wave phenomena and earthquake dynamics. It can be compiled with the __seissol__ make target.
* __MILC__: a collaboration code for lattice QCD calculations, implementing SU(3) lattice gauge theory approaches. It can be compiled with the __milc__ make target.
* __MiniWeather__: a mini-application simulating weather-like flows, reproducing the basic behavior of weather and fluid dynamics workloads. It can be compiled with the __miniw__ make target.
* __GADGET__: a massively parallel code for N-body and hydrodynamical cosmological applications. It can be compiled with the __gadget__ make target.

## External Dependencies

Compiling the HPC applications requires a basic set of generic dependencies. These are __make__, __cmake__, __autoconf__, __pkgconf__ and a functional __LLVM__ and __clang__ toolchain (at least version 15.0, commit 9778ec057cf4, with both source code and binaries). Furthermore, a valid __MPI__ installation is required (OpenMPI is recommended). The specific additional dependencies for each HPC application are compiled automatically as part of the build process, and no user action is required.

Before attempting to compile the HPC applications, make sure to properly configure the __config.mk__ file, as indicated in the main [README](../README.md) document. In order to compile the applications in this section of the reliability benchmark, both a valid OpenBLAS library (__OPENBLAS\_ROOT__ variable) and LLTFI installation (__LLTFI\_ROOT__ variable) are required.

### LLTFI Compilation

The steps to build LLTFI are described in detail on its GitHub [repository](https://github.com/DependableSystemsLab/LLTFI). Please note that a custom version of LLTFI is necessary, including all integration with HDFIT, which can be found in a dedicated __HDFIT__ branch of the repository. For convenience, we also include as part of this repository a git patch and an associated makefile to clone LLTFI, apply all HDFIT-specific changes and then build the toolchain automatically. These files can be found under the __lltfi__ directory. Compilation can be simply achieved as follows:

```
cd lltfi
make
```

Users should point to the appropriate LLVM (at least version 15.0, as mentioned above) source code and install directories by exporting the __LLVM\_SRC\_ROOT__ and __LLVM\_DST\_ROOT__ environment variables respectively. In addition, a functional __Python 3__ installation, the __PyYAML__ module, __Ninja__ and __libprotoc__ are required. For further details, please refer to the LLTFI documentation.

## Compilation Steps

In order to build all HPC applications, the following command can be used:

```
make all
```

Building specific HPC applications can be achieved by using the corresponding make targets. This repository does not contain the complete source code for each application, but rather only the compilation targets necessary to enable HDFIT use: as such, the build process clones each application's source code from the respective repositories, and subsequently compiles them. Any dependencies are compiled as part of this process.

Please note that it may be necessary to use the __LD\_PRELOAD__ environment variable (set to the path to the HDFIT OpenBLAS .so library file, as well as to the path of the LLTFI .so runtime library file) for HPC applications to work correctly. This is to prevent using a different version of the libraries already present on the system, as well as to address the specific linking mechanisms in each HPC application. In case errors are encountered during the build process, please refer to the README documents in each application's directory for more detailed instructions, and to the respective documentation websites.

The supplied HPC applications have been compiled and tested on Ubuntu (18.04, 20.04 and 22.04), as well as on SLES (12SP5). Compatibility with Windows WSL has been verified as well.

__NOTE__: For compatibility reasons, all application targets enforce compliance with the AVX2 instruction set (i.e., Intel Haswell architecture). If you wish to use the application suite on a system that does not support AVX2, edit the Makefile targets accordingly.

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

## GROMACS

The main dependency for GROMACS is FFTW, which is installed automatically by HDFIT. The following two configurations are supplied:

* __in.aminoacids__: Molecular dynamics simulation of amino acid compounds (__topol.tpr__ file);
* __in.flooding2__: Molecular dynamics simulation for A-form Guanylin with endogenous ligands (__topol.tpr__ file).

We include the original GROMACS configuration files from which the __topol.tpr__ inputs are generated, using the __grompp__ utility. Starting a GROMACS simulation can be achieved with the following commands:

```
cd in.flooding2
mpirun -np 1 ../build/bin/gmx mdrun -s topol.tpr -ntmpi 1 -ntomp 1
```

## SeisSol

The SeisSol application comes with a series of dependencies, including HDF5, NetCDF, Libxsmm and several libraries specific to seismology applications.
These are all installed automatically by HDFIT. It should be noted, however, that SeisSol requires __Python 3__ with the __numpy__ package already at
the compilation stage - furthermore, we encountered compilation errors with certain old versions of __cmake__ (i.e., 3.16.3 and older).
SeisSol comes with the following two example configurations:

* __in.tpv12\_13__: Example using a 60-degree dipping normal fault and depth-dependent initial stress conditions (__parameters.par__ file);
* __in.WP2\_LOH1__: Point-source example for earthquake nucleation using the LOH1 model (__parameters.par__ file).

As part of HDFIT, we compile SeisSol in 4th order mode; this can be changed by altering the HDFIT __Makefile__ appropriately.
Starting a SeisSol run can be achieved with the following:

```
cd in.tpv12_13
mkdir out && cd out
OMP_NUM_THREADS=1 mpirun -np 1 ../../build/SeisSol_Release_shsw_4_elastic ../parameters.par
```

## MILC

This version of MILC does not employ any additional dependencies, save from MPI. The following two configurations are supplied:

* __in.rhmc-hisq__: Run using the Rational Hybrid Monte Carlo method with the HISQ action (__su3\_rhmc\_hisq.f2111.1.sample-in__ file);
* __in.rhmd-hisq__: Run using the Rational Hybrid Molecular Dynamics method with the HISQ action (__su3\_rhmd\_hisq.1.sample-in__ file).

In the context of HDFIT, we target the __su3\_rhmc\_hisq__ and __su3\_rhmd\_hisq__ applications, which reflect state-of-the-art use in the lattice
QCD field. Performing a run of the supplied __su3\_rhmc\_hisq__ configuration, for example, can be achieved as in the following:

```
cd in.rhmc-hisq
mkdir out && cd out
mpirun -np 1 ../ks_imp_rhmc/su3_rhmc_hisq ../su3_rhmc_hisq.f2111.1.sample-in
```

## MiniWeather

The main dependency for MiniWeather is the PNetCDF library, which is installed automatically by HDFIT. Additionally, the netCDF4 Python package is
required in order to run experiments. The application comes with two configurations:

* __in.thermal__: Simulates a rising thermal in neutral atmosphere (__mpi\_thermal__ executable);
* __in.mountain__: Simulates propagation of an horizontal wind over a mountain side (__mpi\_mountain__ executable).

It should be noted that MiniWeather configurations are enforced at compile time: therefore, no external configuration files are used, and altering the
configurations requires editing the HDFIT __Makefile__ itself. Performing a run for the rising thermal configuration can be achieved as follows:

```
cd in.thermal
mpirun -np 1 ../cpp/build/mpi_thermal
```

## GADGET

The GADGET application depends on the GSL, HDF5 and FFTW libraries, which are installed automatically by HDFIT. It should be noted that, in order to
run experiments with the application, the Python h5py package is additionally required. The application comes with two configurations:

* __in.G2-galaxy__: collisionless simulation of two disk galaxies merging into each other (__param.txt__ file, requires __galaxy\_littleendian.dat__);
* __in.G2-gassphere__: simulation of gravitational collapse of a gas sphere (__param.txt__ file, requires __gassphere\_littleendian.dat__).

Please note that the extra data files required by the two configurations, which must be placed alongside the respective __param.txt__ files, must be downloaded separately at [this](https://wwwmpa.mpa-garching.mpg.de/gadget4/example_ics.tar) link. Performing simulations can be achieved as follows:

```
cd in.G2-galaxy
mkdir out && cd out
OMP_NUM_THREADS=1 mpirun -np 1 ../../examples/G2-galaxy/Gadget4 ../param.txt
```
