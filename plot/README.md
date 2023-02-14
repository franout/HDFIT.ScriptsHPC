# HDFIT Plotting Scripts

The scripts in this directory can be used to generate plots of various kinds based on CSV files
produced by HDFIT experiments with HPC applications, specifically using the __HDFIT\_runner.sh__ 
script contained under the __test__ directory. In the following are the main definitions used 
as part of our toolchain:

* The *Program Vulnerability Factor* (PVF) and *Architecture Vulnerability Factor* (AVF) terms are used interchangeably. They both represent the percentage of hardware faults (i.e., the percentage of fault injection runs in HDFIT experiments) that lead to application failures. This definition corresponds to that of the *Architecture Vulnerability Factor* (AVF) in *Architecture Design for Soft Errors* by Mukherjee, Shubu.
* The term *User Observable Error* (UOE) is associated with application runs that terminate under visible error conditions (e.g., crashes, unhandled exceptions, failed runtime checks).
* The term *Silent Data Error* (SDE) is used to characterize application runs that terminate successfully, but producing a numerical output that is different from the expected one. We characterize the magnitude of SDEs using the *Normalized Root Mean Squared Error* (NRMSE) as metric.

## CSV Layout Description

The CSV layout of HDFIT experiment files is as follows:

```
rank, op-cnt expected, op-fi, fi-bit, op-cnt got, test pass cnt, test fail cnt, return code, time, nrmse
```

Where each specific field describes the following:

* __rank__: MPI rank targeted by fault injection.
* __op-cnt expected__: total number of operations for MPI rank 0, as measured during a golden run.
* __op-fi__: ID of the operation (for the chosen MPI rank) that is targeted by fault injection.
* __fi-bit__: bit that is targeted by fault injection (e.g., from 0 to 63 for double-precision floating-point). 
* __op-cnt got__: measured number of operations, for MPI rank 0, on this run.
* __test pass cnt__: number of passed tests for this run. In most cases, this is equal to 1 if the application succeeds, and 0 otherwise.
* __test fail cnt__: number of failed tests for this run. In most cases, this is equal to 1 - test pass cnt.
* __return code__: return code of the application for this run.
* __time__: total run time (in milliseconds) for the application in this run.
* __nrmse__: application-dependent nrmse error metric, computed against the output of a golden run.


More fields may be present in the CSV file for certain types of experiments, for example when simulating
RTL hardware models providing additional telemetry. These are as follows:

* __rtl error__: flag describing whether errors in the circuit were raised or not after fault injection.
* __assign uuid__: ID of the instance-agnostic fault injection signal (or UUID) in the circuit chosen for this run.
* __module chain__: chain of RTL module instance IDs describing the specific section of the circuit chosen for fault injection.

## Available Scripts

The following plotting scripts can be used to plot relevant data for any kind of HDFIT experiment:

* __HDFIT\_plot\_error\_curve.py__: plots the application's failure rate (i.e., PVF or AVF) in function of a NRMSE threshold, as computed from each run's output. It additionally prints a variety of failure statistics.
* __HDFIT\_plot\_per\_bit\_bar.py__: produces a bar plot of the application's failure rate (i.e., PVF or AVF) in function of injection bit positions.
* __HDFIT\_plot\_heatmap.py__: produces a heatmap of the failure rate (i.e., PVF or AVF) in function of both injection bit positions, and code positions throughout the application's runtime.
* __HDFIT\_plot\_times.py__: produces a bar plot of the percentage overhead, in function of injection bit positions, with respect to the reference number of operations (as computed by HDFIT).

The following scripts, on the other hand, plot data that is specific to experiments performed using RTL-level (i.e., Netlist) simulations:

* __HDFIT\_plot\_uuid.py__: produces a bar plot of failure occurrences in function of UUIDs (i.e., the IDs of the circuit signals to which faults are injected).
* __HDFIT\_plot\_modules.py__: produces a bar plot of failure occurrences in function of RTL module names. This script can additionally augment the input CSV file with RTL module name information.
* __HDFIT\_plot\_instance\_tree.py__: produces a breakdown of failure occurrences according to the actual hierarchy of instantiated RTL components. The resulting instance tree can be visualized in text form or plotted using Graphviz.

## Overall Usage

At a minimum, all scripts need as input argument the filename of the CSV file to be processed:

```
python3 HDFIT_plot_error_curve.py ../../cp2k/out.C2H4/HDFIT-CP2K-C2H4-29.08.2022-transient.csv
```

The __HDFIT\_plot\_modules.py__ and __HDFIT\_plot\_instance\_tree.py__ scripts additionally require as input the path to the C++ source file which contains 
the HDFIT specification for fault injection signals. Additional configuration parameters are available, which can be discovered using each script's 
`-h` option. All scripts are meant to be used in a __Python 3__ environment and require the __numpy__, __matplotlib__ and __seaborn__ packages.
