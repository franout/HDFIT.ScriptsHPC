#!/bin/bash
# Copyright (C) 2022 Intel Corporation
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License, as published
# by the Free Software Foundation; either version 3 of the License,
# or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# ----- ENVIRONMENT SETUP -----
echo "------ HPC application runner for HDFIT ------"
# Populating environment variables from the test script
source $1
# Exporting proper output redirection environment variable
export BLASFI_OUTPUT=$FI_STREAM
# Base directory of the application
export FI_BASEDIR=$2

if $3 ; then 
golden_simulation=true
else 
golden_simulation=false
fi 

# Output directory 
export FI_RESDIR_SUP=$FI_BASEDIR/"out."$FI_CONFNAME
# Input directory
export FI_CONFDIR=$FI_BASEDIR/"in."$FI_CONFNAME
# Experiment name
export FI_TASKNAME="HDFIT_${FI_APPNAME}_${FI_CONFNAME}_$(date +%d-%m-%Y)"
# Parsing the paths to required libraries from the parent config.mk file (or check env variables)
if [ -n "${OPENBLAS_PATH}" ] ; then 
OPENBLAS_PATH=$(cat $FI_THISDIR/../config.mk | grep -aoP "(?<=OPENBLAS_ROOT = ).*")
else
OPENBLAS_PATH=${OPENBLAS_ROOT}
fi 

if [ -n "${LLTFI_PATH}" ] ; then 
LLTFI_PATH=$(cat $FI_THISDIR/../config.mk | grep -aoP "(?<=LLTFI_ROOT = ).*")
else
LLTFI_PATH=${LLTFI_ROOT}
fi 
export FI_PRELOAD="$OPENBLAS_PATH/libopenblas.so:$LLTFI_PATH/runtime_lib/libllfi-rt.so"

# Performs special pre-processing for some applications' input
doSpecialInput()
{
	# HPL expects an HPL.dat file to be present in the cwd - so we need to
	# create a symlink to said file in the upper directory
	if [ $FI_APPNAME = "HPL" ]; then
		ln -s ${inputs_dir}/in.${FI_CONFNAME}/HPL.dat HPL.dat 2>/dev/null
	fi
}

# Performs special post-processing for some applications' output
doSpecialOutput()
{
        # The hack of all hacks to print DFT energy information in a sensible manner
        # Coded ad-hoc for the NWCHEM 3carbo_dft input with 16 atoms
        if [ $FI_APPNAME = "NWCHEM" ] && [ $FI_CONFNAME = "3carbo" ]; then
                cat $FI_RESDIR/run$1.log | grep -a -A20 "DFT ENERGY GRADIENTS" > $FI_RESDIR/Bzisox_qmd.dft
        # MILC only writes to stdout - creating an ad-hoc output file
        elif [ $FI_APPNAME = "MILC" ]; then
                cat $FI_RESDIR/run$1.log > $FI_RESDIR/$FI_CONFNAME.sample-out
        # Converting the GROMACS binary format to text
        elif [ $FI_APPNAME = "GROMACS" ]; then
                LD_PRELOAD=$FI_PRELOAD $FI_COMMAND dump -f traj.trr > traj.txt 2>/dev/null
        fi
}

# Performs a single FI (or golden) application run
doApplicationRun()
{
	cd $FI_CONFDIR/run$2
	mkdir -p -m u=rwX,go=rX $FI_RESDIR/run$1.DAT
	doSpecialInput $1

	# Determining range of CPU cores to assign to this run
	if [ $FI_CPUMAP -gt 0 ] && [ $FI_PARRUNS -gt 1 ]
	then
		CPU_FIRST=$(( $2 * $FI_MPIRANKS ))
		CPU_LAST=$(( ($2 + 1) * $FI_MPIRANKS - 1 ))
		CPU_MAP="-cpu-list $(seq -s, $CPU_FIRST 1 $CPU_LAST)"
	else
		CPU_MAP="--map-by core"
	fi
	export evcd_file_path=${SDENV_VULNERABILITY_INPUT_DIR}"/run$1_${FI_APPNAME}_${FI_CONFNAME}.evcd"
	TIMESTART=$(date +%s%3N)
	if ${golden_simulation} ; then 
	## no timeout 
	echo "mpirun -x LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$FI_THISDIR/../apps/deps/install/lib \
		       -x LD_PRELOAD=$FI_PRELOAD -np $FI_MPIRANKS $CPU_MAP $FI_COMMAND $FI_INPUT > $FI_RESDIR/run$1.log 2>&1"
	mpirun -x LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$FI_THISDIR/../apps/deps/install/lib \
		       -x LD_PRELOAD=$FI_PRELOAD -np $FI_MPIRANKS $CPU_MAP $FI_COMMAND $FI_INPUT > $FI_RESDIR/run$1.log 2>&1
	else 
	echo "timeout -k 60s $FI_TIMEOUT mpirun -x LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$FI_THISDIR/../apps/deps/install/lib \
		       -x LD_PRELOAD=$FI_PRELOAD -np $FI_MPIRANKS $CPU_MAP $FI_COMMAND $FI_INPUT > $FI_RESDIR/run$1.log"
	timeout -k 60s $FI_TIMEOUT mpirun -x LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$FI_THISDIR/../apps/deps/install/lib \
		       -x LD_PRELOAD=$FI_PRELOAD -np $FI_MPIRANKS $CPU_MAP $FI_COMMAND $FI_INPUT > $FI_RESDIR/run$1.log 2>&1
	fi 
	FI_RETCODE=$?
	echo $FI_RETCODE > $FI_RESDIR/run$1.retcode
	echo $(( $(date +%s%3N) - $TIMESTART )) > $FI_RESDIR/run$1.time
	cat $FI_RESDIR/run$1.log | grep -a "\[HDFIT\]" > $FI_RESDIR/run$1.hdfit
	if [ ! -z "$FI_OUTPUT" ]; then
		doSpecialOutput $1
		mv $FI_OUTPUT $FI_RESDIR/run$1.DAT/ 2>/dev/null
	fi
	# Wiping out directory contents just in case
	rm -rf $FI_CONFDIR/run$2/*

	# Performing logging
	if [[ $FI_RETCODE -ne 0 ]]
	then
		echo "    Run no. $1 failed with error code $FI_RETCODE. " \
			  "Check $FI_RESDIR/run$1.log for details."
	else
		echo "    Run no. $1 completed successfully."
	fi
}

# Experiment runner corresponding to a single process
doExperimentJob()
{
	JOB_NUMRUNS=$1
	JOB_PARRUNS=$2
	JOB_BASEIDX=$3
	JOB_INDEX=$4
	# Introducing noise before starting experiment
	sleep 0.$(( $JOB_INDEX + 1 ))
        for (( i=$JOB_BASEIDX; i<$JOB_NUMRUNS; i++ ))
        do
		REALRUNID=$(( $i*$JOB_PARRUNS + $JOB_INDEX))
		doApplicationRun $REALRUNID $JOB_INDEX
        done
}

# Overall experiment runner wrapper
doExperiment()
{
	EXP_NUMRUNS=$1
	EXP_PARRUNS=$2
	EXP_BASEIDX=0
	mkdir -p -m u=rwX,go=rX $FI_RESDIR
	if [[ $FI_RECOVER -gt 0 ]]
	then
		EXP_BASEIDX=$(( $(ls -l $FI_RESDIR/*.retcode 2>/dev/null | wc -l) / $EXP_PARRUNS ))
	fi

	# Trapping signals for (somewhat) clean termination
	trap "trap - TERM && kill -- -$$ 2>/dev/null" INT TERM
	# Calling user-defined experiment runner
	doExperimentLauncher $EXP_NUMRUNS $EXP_PARRUNS $EXP_BASEIDX
}

# ----- INITIALIZATION -----
# Cleaning up old results
if [[ $FI_RECOVER -eq 0 ]]
then
	rm -rf $FI_RESDIR_SUP
fi
mkdir -p -m u=rwX,go=rX $FI_RESDIR_SUP

for (( j=0; j<$FI_PARRUNS; j++))
do
        mkdir -p -m u=rwX,go=rX $FI_CONFDIR/"run"$j
done

echo "    Application: $FI_APPNAME"
echo "    Configuration: $FI_CONFNAME"
echo "    Number of iterations: $FI_NUMRUNS"
echo "    Number of runs per iteration: $FI_PARRUNS"
echo "    Number of MPI processes: $FI_MPIRANKS"
echo "    Number of OMP threads: $OMP_NUM_THREADS"
echo "    Run timeout value: $FI_TIMEOUT"
echo "    Output stream: $FI_STREAM"
echo "    Experiment location: $FI_RESDIR_SUP"
echo "----------------------------------------------"

if ${golden_simulation}; then 
# ----- GOLDEN RUN -----
export BLASFI_MODE="NONE"
export BLASFI_CORRUPTION="FLIP"
export BLASFI_BITS="EVERYWHERE"
export BLASFI_OPSCNT=1
export FI_RESDIR=$FI_RESDIR_SUP/"golden"

echo "Performing golden run..."
export evcd_file_path=$FI_RESDIR"/golden_run.evcd"

doExperiment 1  1 # $FI_NUMRUNS $FI_PARRUNS

# Checking outcome of golden run
GOLDEN_RETCODE=$(cat $FI_RESDIR/run0.retcode 2>/dev/null)
if [[ ! -e $FI_RESDIR/run0.retcode ]]; then
	echo "    ERROR: The golden run could not be executed. Cannot continue."
	exit 0
elif [[ $GOLDEN_RETCODE -ne 0 ]]; then
	echo "    ERROR: The golden run failed with error code $GOLDEN_RETCODE. Cannot continue."
	echo "    Check $FI_RESDIR/run0.log for details."
	exit 0
fi

fi 

## TO BE IMPLEMENTED 
if false ; then 
doExperiment 1 1 

# Checking outcome of golden run
GOLDEN_RETCODE=$(cat $FI_RESDIR/run0.retcode 2>/dev/null)
if [[ ! -e $FI_RESDIR/run0.retcode ]]; then
	echo "    ERROR: The golden run could not be executed. Cannot continue."
	exit 0
elif [[ $GOLDEN_RETCODE -ne 0 ]]; then
	echo "    ERROR: The golden run failed with error code $GOLDEN_RETCODE. Cannot continue."
	echo "    Check $FI_RESDIR/run0.log for details."
	exit 0
fi

# Extracting number of BLAS operations from golden run
export BLASFI_OPSCNT=$(cat $FI_RESDIR/run0.log | grep -aoP "(?<=Rank 0: OpsCnt = )[0-9]+")
# Resolving wildcards and computing expanded list of output files
if [ ! -z "$FI_OUTPUT" ]; then
	export FI_OUTPUT_EXP=$(cd $FI_RESDIR/run0.DAT && ls -d $FI_OUTPUT | tr "\n" " ")
fi


# ----- TRANSIENT FAULTS -----
export BLASFI_MODE="TRANSIENT"
export BLASFI_CORRUPTION="FLIP"
export BLASFI_BITS="EVERYWHERE"
export FI_RESDIR="$FI_RESDIR_SUP/fi-transient"

echo "Performing transient faults experiment..."
doExperiment $FI_NUMRUNS $FI_PARRUNS
echo "Computing summary CSV file..."
echo "${PYTHON} $FI_THISDIR/HDFIT_computeCSV.py $FI_APPNAME $FI_RESDIR $BLASFI_OPSCNT $FI_OUTPUT_EXP > $FI_RESDIR_SUP/$FI_TASKNAME.csv"
${PYTHON} $FI_THISDIR/HDFIT_computeCSV.py $FI_APPNAME $FI_RESDIR $BLASFI_OPSCNT $FI_OUTPUT_EXP > "$FI_RESDIR_SUP/$FI_TASKNAME.csv"

fi 