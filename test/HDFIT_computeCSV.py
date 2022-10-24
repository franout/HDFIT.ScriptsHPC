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

from HDFIT_parsers import *
import os, sys
import numpy as np

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Insufficient number of arguments!")
        print("Usage: HDFIT-computeCSV.py <APPNAME> <RESDIR> <OPSCNT> [<OUTPUT 0> <OUTPUT ...> <OUTPUT N>]")
        exit(1)

    appName = sys.argv[1]
    homeDir = sys.argv[2]
    opsCntE = sys.argv[3]
    fiDir   = homeDir
    golDir  = homeDir + "/../golden/"
    outputAvail = len(sys.argv) > 4

    # In this section we compute the max-min bounds of the data for the NRMSE computation
    # The formula we use is NRMSE = sqrt(MSE) / (max - min)
    # Normalization is performed against the max-min range in order to represent the
    # error perceived by HPC users when consuming application output data
    if outputAvail:
        goMats  = loadDmp([golDir + "run0.DAT/" + d for d in sys.argv[4:]], appName)
        goMax   = [ np.max(mat) for mat in goMats ]
        goMin   = [ np.min(mat) for mat in goMats ]

        goConstant = False
        for k in range(len(goMin)):
            if goMax[k] == goMin[k]:
                goMin[k] = 0
                goMax[k] = np.abs(goMax[k]) if goMax[k]!=0 else 1
                goConstant = True

        if goConstant:
            print("    One or more parts of the golden output are constant.", file=sys.stderr)
            print("    Using max as NRMSE normalization factor.", file=sys.stderr)

    # Identifying available data for this experiment
    headerStr   = "rank, op-cnt expected, op-fi, fi-bit, op-cnt got, test pass cnt, test fail cnt, return code, time, nrmse"
    containsRTL = probeRTLInfo(fiDir + "/run0.hdfit") 
    if containsRTL:
        headerStr = headerStr + ", rtl error, assign uuid, module chain"
    print(headerStr)

    numRuns = len([el for el in os.listdir(fiDir) if ".retcode" in el]) 
    for i in range(numRuns):
        namePrefix = fiDir + "/run"
        try:
            firetcodeF = open(namePrefix + str(i) + ".retcode")
            fitimeF    = open(namePrefix + str(i) + ".time")
            if outputAvail:
                fiMats = loadDmp([namePrefix + str(i) + ".DAT/" + d for d in sys.argv[4:]], appName)
            rank, fiOp, fiBit, opsCnt = loadFaultInfo(namePrefix + str(i) + ".hdfit")
            if containsRTL:
                rtlError, assignUUID, moduleChain = loadRTLInfo(namePrefix + str(i) + ".hdfit")
        except IOError:
            print("    Run %d terminated unexpectedly!" % (i), file=sys.stderr)
            continue

        firetcodeL = firetcodeF.readline().strip()
        fitimeL = fitimeF.readline().strip()
        if not firetcodeL or not fitimeL:
            print("    Run %d has missing fault information!" % (i), file=sys.stderr)
            continue
        
        if not outputAvail:
            diffCoeff = 0.0
        elif len(goMats) != len(fiMats) or any(fiMats[idx].shape != goMats[idx].shape for idx in range(len(goMats))):
            print("    Unexpected output shape for run %d compared to golden run." % (i), file=sys.stderr)
            diffCoeff = np.inf
        else:
            try:
                diffCoeffs = [rmse(goMats[idx], fiMats[idx]) / (goMax[idx] - goMin[idx]) for idx in range(len(goMats))]
                diffCoeff  = np.average(diffCoeffs) * 100
            except ValueError:
                print("    Error while computing the NRMSE for run %d!" % (i), file=sys.stderr)
                diffCoeff = np.inf
        
        diffCoeff  = str(diffCoeff)
        passcnt    = "1" if firetcodeL=="0" else "0"
        failcnt    = str(1 - int(passcnt))
        outStr     = ",".join([rank, opsCntE, fiOp, fiBit, opsCnt, passcnt, failcnt, firetcodeL, fitimeL, diffCoeff])
        if containsRTL:
            outStr = outStr + "," + ",".join([rtlError, assignUUID, moduleChain])
        print(outStr)

    exit(0)
