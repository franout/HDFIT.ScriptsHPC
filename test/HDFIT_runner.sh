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

# Checking arguments
if [ "$#" -ne 1 ]
then
  echo "Usage: $0 <TEST SCRIPT>" >&2
  exit 1
fi

# Base directory of this script
export FI_THISDIR=$(dirname $(realpath "$BASH_SOURCE"))
# Full path to OpenBLAS .so file - necessary for some applications
export FI_PRELOAD=PATH_TO_CUSTOM_OPENBLAS_LIB/libopenblas.so
# Activates experiment recovery (e.g., if an HPC job fails)
export FI_RECOVER=0
# Enables automatic CPU core pinning (works only for single-node jobs)
export FI_CPUMAP=1

# Simple user-defined launcher for multiple experiment jobs
doExperimentLauncher()
{
	RUN_NUMRUNS=$1
	RUN_PARRUNS=$2
	RUN_BASEIDX=$3
        for (( j=0; j<$RUN_PARRUNS; j++))
        do
       		doExperimentJob $RUN_NUMRUNS $RUN_PARRUNS $RUN_BASEIDX $j &
	done
	echo "    Spawned $RUN_PARRUNS parallel runners. Waiting for completion..."
	wait
	echo "    All parallel experiment runners terminated."
}

# Everything ready - launch internal runner script
source $FI_THISDIR/HDFIT_runner_internal.sh "$@"
