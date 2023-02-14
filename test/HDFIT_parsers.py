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

import os, sys
import re
import numpy as np
import xml.etree.ElementTree as et


def loadDmp(fPaths, appName):
    """
    Generic wrapper for loading application data.
  
    This function is a generic wrapper for loading data associated with
    HDFIT-supported applications. Data is always returned as a list of
    matrices, one for each metric made available.
  
    Parameters:
        fPaths (list[string]): Paths to files to be parsed
  
    Returns:
        list[np.array]: Matrices with parsed data
    """
    try:
        if appName=="QE":
            return loadDmpQE(fPaths)
        elif appName=="NWCHEM":
            return loadDmpNWCHEM(fPaths)
        elif appName=="REMHOS":
            return loadDmpREMHOS(fPaths)
        elif appName=="CP2K":
            return loadDmpCP2K(fPaths)
        elif appName=="QMCPACK":
            return loadDmpQMCPACK(fPaths)
        elif appName=="LAMMPS":
            return loadDmpLAMMPS(fPaths)
        elif appName=="MILC":
            return loadDmpMILC(fPaths)
        elif appName=="MINIW":
            return loadDmpMINIW(fPaths)
        elif appName=="SEISSOL":
            return loadDmpSEISSOL(fPaths)
        elif appName=="GROMACS":
            return loadDmpGROMACS(fPaths)
        elif appName=="GADGET":
            return loadDmpGADGET(fPaths)
        else:
            raise ValueError("Application type unknown!")
    except IOError:
        return [np.zeros((1, 1))]

def rmse(y_true, y):
    """
    Function to compute the RMSE.
  
    Parameters:
        y_true (np.array): Numpy array of true values
        y (np.array): Numpy array of estimated values
  
    Returns:
        float: Root mean square error value
    """
    return np.linalg.norm(y - y_true) / np.sqrt(y_true.size)

def parseFloats(splitLine):
    """
    Converts a list of strings into a list of floats.
  
    Parameters:
        splitLine (list[string]): Strings to be converted
  
    Returns:
        list[float]: Converted values
    """
    tmpVec = []
    for el in splitLine:
        try:
            tmpVec.append(float(el))
        except ValueError:
            print("    Cannot parse floating-point value %s. Treating as 0." % (el), file=sys.stderr)
            tmpVec.append(0)
    return tmpVec

def loadFaultInfo(fPath):
    """
    Parses basic HDFIT metrics from a log file.
    
    This function will look for and parse basic HDFIT metrics
    within the chosen file. The metrics include the ops count,
    as well as the MPI rank, op and bit of fault injection. These
    are available for all kinds of HDFIT experiments.    
  
    Parameters:
        fPath (string): Path to file to be parsed
  
    Returns:
        string: MPI rank of fault injection
        string: Operation ID of fault injection
        string: Bit of fault injection
        string: Total operations counter
    """
    f = open(fPath, "r")
    # This will throw an IOError if fPath does not exist
    fText = f.read()
    f.close()

    rank   = re.search("(?<=FI enabled on rank = )[0-9]+(?=\n)", fText)
    fiop   = re.search("(?<=FI at op = )[0-9]+(?=\n)", fText)
    fiBit  = re.search("(?<=Bit pos = )[0-9]+(?=\n)", fText)

    if rank is None or fiop is None or fiBit is None:
        raise IOError("Missing fault information in file %s!" % fPath)

    opsCnt0 = re.search("(?<=Rank 0: OpsCnt = )[0-9]+(?=\n)", fText)
    opsCntR = re.search("(?<=Rank " + rank.group(0) + ": OpsCnt = )[0-9]+(?=\n)", fText)
    # Picking as total ops count, in this order:
    # 1) Ops count of rank 0, if available
    # 2) Ops count of FI rank, if available
    # 3) Ops count at time of FI, on FI rank
    opsCnt = opsCnt0 if opsCnt0 is not None else opsCntR if opsCntR is not None else fiop

    if opsCnt is None:
        raise IOError("Missing ops count information in file %s!" % fPath)

    return rank.group(0), fiop.group(0), fiBit.group(0), opsCnt.group(0)

def loadRTLInfo(fPath):
    """
    Parses RTL-specific HDFIT metrics from a log file.
    
    This function will parse HDFIT metrics that are specific
    to experiments performed using RTL or Netlist-level hardware
    simulations. These include the signal ID of fault injection, the
    associated module instance chain, and a flag to signal RTL errors.
  
    Parameters:
        fPath (string): Path to file to be parsed
  
    Returns:
        string: RTL error status
        string: Signal ID of fault injection
        string: Module instance chain of fault injection
    """
    f = open(fPath, "r")
    fText = f.read()

    rtlError    = re.search("(?<=RTL errors = )[0-9]+(?=\n)", fText)
    assignUUID  = re.search("(?<=Assign UUID = )[0-9]+(?=\n)", fText)
    moduleChain = re.search("(?<=Module instance chain = ).+(?=\n)", fText)
    f.close()

    if rtlError is None or assignUUID is None or moduleChain is None:
        raise IOError("Missing RTL information in file %s!" % fPath)

    return rtlError.group(0), assignUUID.group(0), moduleChain.group(0)

def probeRTLInfo(fPath):
    """
    Determines the presence of RTL data in a log file.
    
    This utility function simply determines whether a HDFIT log
    file contains RTL-specific data or not. It does not return the
    data itself.
  
    Parameters:
        fPath (string): Path to file to be probed
  
    Returns:
        bool: True if RTL data is present, False otherwise
    """
    try:
        loadRTLInfo(fPath)
    except IOError:
        return False
    return True

# ------------------------------------ HPL Parser ------------------------------------

def loadResidualsHPL(fPath):
    """
    A parser for HPL log data.

    This function extracts the residual value produced within a single HPL benchmark
    log - unlike the other parsers, it is not able to process multiple logs at once.
    The function is not used by default in the HDFIT_computeCSV.py script, but can
    be easily integrated to add HPL residual information to the output CSV file.

    Parameters:
        fPath (string): Path to HPL log file

    Returns:
        string: Parsed residual value in string form
    """
    f = open(fPath, "r")
    fText = f.read()
    f.close()

    resLine   = re.search("\*N\)\=.+(?=\n)", fText)

    if resLine is None:
        return "inf"

    resVal  = re.search("(?<==   )\s*[0-9ena\.\+\-]+", resLine.group(0))
    return resVal.group(0).strip()

# ---------------------------------- GADGET Parser -----------------------------------

def loadDmpGADGET(fPaths):
    """
    A parser for GADGET output data.

    This function uses the h5py module to parse the output of
    HDF5 files produced by GADGET simulations. This consists
    of multiple metrics, whose IDs are retrieved dynamically.

    Parameters:
        fPaths (list[string]): Paths to files to be parsed

    Returns:
        list[np.array]: Matrices with parsed data
    """
    # h5py required to parse output files
    import h5py
    dmpMats = {}

    for fp in fPaths:
        f = h5py.File(fp, "r")

        # Parsing and concatenating all specified netCDF files
        for part in f.keys():
            if "PartType" in part:
                for met in f[part].keys():
                    if "ParticleIDs" not in met:
                        metID = part + met
                        tmpMet = np.reshape(f[part][met][:].flatten(), (1, -1))
                        dmpMats[metID] = np.copy(tmpMet) if metID not in dmpMats else np.concatenate((dmpMats[metID], tmpMet), axis=0)

    return [dmpMats[met] if met in dmpMats else np.zeros((1,1)) for met in dmpMats.keys()]

# ---------------------------------- GROMACS Parser ----------------------------------

def loadDmpGROMACS(fPaths):
    """
    A parser for GROMACS output data.

    This function currently supports GROMACS trajectory data that has been converted
    to text form using the gmx dump utility. As such, position, velocity and force
    data is parsed - energy data is not currently supported.

    Parameters:
        fPaths (list[string]): Paths to files to be parsed

    Returns:
        list[np.array]: Matrices with parsed data
    """
# ---------- Helper function
    def parseLineGROMACS(line):
        line = line.split('{')[1][0:-1]
        splitLine = [el.strip() for el in line.split(",") if el.strip()!=""]
        return parseFloats(splitLine)
# ----------

    metrics = ['x[', 'v[', 'f[']
    dmpMats = {}
    for met in metrics:
        dmpMats[met] = []

    for fp in fPaths:
        if ".txt" in fp:
            f = open(fp, 'r')
            for line in f:
                line = line.strip()
                met = line[0:2]
                if met in metrics:
                    dmpMats[met].append(parseLineGROMACS(line))

    return [np.array(dmpMats[met], ndmin=2) if len(dmpMats[met]) > 0 else np.zeros((1,1)) for met in metrics]

# ------------------------------------ SEISSOL Parser ------------------------------------

def loadDmpSEISSOL(fPaths):
    """
    A parser for SeisSol data.

    At the moment, this function parses the content of SeisSol
    receiver files. When multiple files are specified, their
    contents are simply concatenated, differentiating between
    normal "receiver" files and "faultreceiver" ones. This leads
    to one large matrix containing data from all receivers: this
    is then split column-wise, returning each metric separately.

    Parameters:
        fPaths (list[string]): Paths to files to be parsed

    Returns:
        list[np.array]: Matrices with parsed data
    """
# ---------- Helper function
    def parseReceiverSEISSOL(fPathInt):
        f = open(fPathInt, 'r')
        tmpMat = []
        lineCnt = 0

        for line in f:
            lineCnt = lineCnt + 1
            if lineCnt > 2 and line[0] != "#":
                splitLine = [el.strip() for el in line.split(" ") if el.strip()!=""]
                # Dropping the timestep field
                tmpMat.append(parseFloats(splitLine[1:]))

        if len(tmpMat)==0:
            raise IOError("Empty SEISSOL receiver file!")

        # Discarding index column
        return np.array(tmpMat, ndmin=2)

    def parseXDMFSEISSOL(fPathInt):
        # Note: the XDMF parser for SEISSOL requires an ad-hoc Python module
        # and is therefore disabled by default. Including XDMF surface files
        # did not yield significant NRMSE differences during testing
        import seissolxdmf
        surDat = {}
        sx = seissolxdmf.seissolxdmf(fPathInt)
        for met in sx.ReadAvailableDataFields():
            if met != "partition":
                surDat[met] = np.reshape(sx.ReadData(met).flatten(), (1, -1))
        return surDat
# ----------
    recDat = None
    fauDat = None
    surDat = []

    for fp in fPaths:
        if "faultreceiver" in fp:
            recTmp = parseReceiverSEISSOL(fp)
            fauDat = np.copy(recTmp) if fauDat is None else np.concatenate((fauDat, recTmp), axis=0)
        elif "receiver" in fp:
            recTmp = parseReceiverSEISSOL(fp)
            recDat = np.copy(recTmp) if recDat is None else np.concatenate((recDat, recTmp), axis=0)
        # Below branch disabled by default
        elif "xdmf" in fp and False:
            surTmp = parseXDMFSEISSOL(fp)
            surDat.append(surTmp)

    # Separating single columns of the final matrix
    outMat = [recDat[:, idx] for idx in range(recDat.shape[1])]
    if fauDat is not None:
        outMat = outMat + [fauDat[:, idx] for idx in range(fauDat.shape[1])]
    if len(surDat) > 0:
        outMat = outMat + [tmpDat[met] for tmpDat in surDat for met in tmpDat.keys()]

    return outMat

# ---------------------------------- MINIW Parser -----------------------------------

def loadDmpMINIW(fPaths):
    """
    A parser for MiniWeather output data.

    This function uses the netCDF4 module to parse the output of
    MiniWeather runs. This consists of several metrics (indicated by
    the metrics list), each sampled over a t,x,z domain. The x,z images
    are concatenated over all timestamps (as this does not impact NRMSE
    computation anyway) for each metric separately, and then returned.

    Parameters:
        fPaths (list[string]): Paths to files to be parsed

    Returns:
        list[np.array]: Matrices with parsed data
    """
    # netCDF required to parse MiniWeather output files
    import netCDF4
    dmpMats = {}
    metrics = ['dens', 'uwnd', 'wwnd', 'theta']

    for fp in fPaths:
        nc = netCDF4.Dataset(fp)

        # Parsing and concatenating all specified netCDF files
        for met in metrics:
            tmpMet = np.reshape(nc.variables[met][:].flatten(), (1, -1))
            dmpMats[met] = np.copy(tmpMet) if met not in dmpMats else np.concatenate((dmpMats[met], tmpMet), axis=0)

    return [dmpMats[met] if met in dmpMats else np.zeros((1,1)) for met in metrics]

# ----------------------------------- MILC Parser -----------------------------------

def loadDmpMILC(fPaths):
    """
    A parser for MILC output data.

    This function parses relevant metrics generated by MILC's su3_rhmc_hisq and
    su3_rhmd_hisq executables. The specific list of metrics is indicated by the
    dataDict and idxDict dictionaries - after parsing, each metric ends up in
    a separate matrix.

    Parameters:
        fPaths (list[string]): Paths to files to be parsed

    Returns:
        list[np.array]: Matrices with parsed data
    """

# ---------- Helper function
    def parseLineMILC(line, idx):
        line = line.replace("\t", " ")
        splitLine = [el.strip() for el in line.split(" ") if el.strip()!=""]
        return parseFloats(splitLine[idx[0]:idx[1]])
# ----------
    dmpMats = []
    # Dictionary of supported metrics for parsing
    dataDict = {"G_LOOP" : [], "P_LOOP" : [], "ACTION" : [], "GACTION" : [], 
                "FACTION" : [], "PBP" : [], "PLAQ" : []}
    # Dictionary of indices to use for parsing lines
    idxDict = {"G_LOOP" : (4, 5), "P_LOOP" : (1, 3), "ACTION" : (3, 7), "GACTION" : (1, 2),
               "FACTION" : (4, 5), "PBP" : (3,7), "PLAQ" : (1, 3)}

    for fp in fPaths:
        f = open(fp, 'r')
        for line in f:
            tok = line.split(":")[0]
            if tok in dataDict:
                dataDict[tok].append(parseLineMILC(line, idxDict[tok]))

    for met in dataDict:
        dmpMats.append(np.array(dataDict[met], ndmin=2) if len(dataDict[met])>0 else np.zeros((1,1)))
    return dmpMats

# ------------------------------------ QE Parser ------------------------------------

def loadDmpQE(fPaths):
    """
    A parser for Quantum Espresso data.
    
    This function is able to parse QE data produced by the PW application,
    in xml format. In particular, it mainly parses the ks_energies values
    (plus a couple other metrics), which are commonly found in SCF runs.
    Each metric (e.g., ks_energies and forces) is stored in a separate matrix.
  
    Parameters:
        fPaths (list[string]): Paths to files to be parsed
  
    Returns:
        list[np.array]: Matrices with parsed data
    """
# ---------- Helper function
    def parseLineQE(line):
        line = line.replace("\n", " ")
        splitLine = [el.strip() for el in line.split(" ") if el.strip()!=""]
        return parseFloats(splitLine)
# ----------
    dmpMats = []

    for fp in fPaths:
        tree = et.parse(fp)
        root = tree.getroot()

        dmpEig = []
        for ks in root.iter('ks_energies'):
            dmpEig.append(parseLineQE(ks.find("eigenvalues").text))
        dmpMats.append(np.array(dmpEig, ndmin=2))

        dmpOcc = []
        for ks in root.iter('ks_energies'): 
            dmpOcc.append(parseLineQE(ks.find("occupations").text))
        dmpMats.append(np.array(dmpOcc, ndmin=2))

        fs = root.find('output/forces')
        if fs is not None:
            dmpFs = [parseLineQE(fs.text)]
            dmpMats.append(np.array(dmpFs, ndmin=2))

        hs = root.find('output/dft/dftU/Hubbard_ns_mod')
        if hs is not None:
            dmpHs = [parseLineQE(hs.text)]
            dmpMats.append(np.array(dmpHs, ndmin=2))

    return dmpMats

# ------------------------------------ NWCHEM Parser ------------------------------------

def loadDmpNWCHEM(fPaths):
    """
    A parser for NWChem data.
    
    This function parses output data produced by NWChem. In particular,
    it is able to parse .emotion and .xyz files produced as part of
    Car-Parrinello simulations, plus .trj files that are produced from
    MD simulations. In addition, the function is able to parse DFT energy
    gradient values - this, however, refers to a custom output format that 
    is used in HDFIT experiments; normally, DFT energy gradients are only 
    found as part of the main NWChem log. Each metric type is stored into 
    a separate matrix. It should be noted that coordinates and gradients 
    in XYZ and DFT files are split into two separate matrices.
  
    Parameters:
        fPaths (list[string]): Paths to files to be parsed
  
    Returns:
        list[np.array]: Matrices with parsed data
    """
    dmpMats = []

    for fp in fPaths:
        f = open(fp, 'r')
        isXYZ = ".xyz" in fp
        isDFT = ".dft" in fp
        tmpMat = []
        
        if isXYZ or isDFT or ".emotion" in fp:
            for line in f:
                splitLine = [el.strip() for el in line.split(" ") if el!=""]
                # Skipping empty or header lines
                # DFT lines have 8 fields, 6 of which numerical
                if (not isDFT and len(splitLine) > 1) or (isDFT and len(splitLine) == 8):
                    splitLine = splitLine[1:] if not isDFT else splitLine[2:]
                    tmpMat.append(parseFloats(splitLine))

            # Splitting ion positions and velocities in two separate matrices
            # if we are reading a .xyz or .dft file
            tmpMat = np.array(tmpMat, ndmin=2)
            dmpMats = dmpMats + ([tmpMat[:,0:3], tmpMat[:,3:6]] if isXYZ or isDFT else [tmpMat])
        elif ".trj" in fp:
            validLines = False
            for line in f:
                splitLine = [el.strip() for el in line.split(" ") if el!=""]
                if splitLine[0] == "TFTF":
                    validLines = True
                elif splitLine[0] == "frame":
                    validLines = False
                elif validLines:
                    tmpMat.append(parseFloats(splitLine))

            dmpMats.append(np.array(tmpMat, ndmin=2))
        else:
            raise IOError("Unknown file extension for parsing!")

    return dmpMats

# ------------------------------------ REMHOS Parser ------------------------------------

def loadDmpREMHOS(fPaths):
    """
    A parser for Remhos data.
    
    This function parses MFEM output data produced by Remhos. In particular,
    it is able to parse both mesh and solution files across multiple time steps.
    It should be noted that the parser expects output files that use a specific
    type of ordering (type 0), which is the one used across our test cases, and
    may not work with other types. The parser returns two matrices, one for mesh
    and the other for solution data, with rows within each matrix corresponding to
    different time steps.
  
    Parameters:
        fPaths (list[string]): Paths to time step directories containing MFEM data
  
    Returns:
        list[np.array]: Matrices with parsed data
    """
# ---------- Helper function
    def parseFileREMHOS(fPathInt):
        f = open(fPathInt, 'r')
        dataVec = []
        ordFound = False
        prevLine = ""

        for line in f:
            if prevLine.strip() == "Ordering: 0":
                ordFound=True
            elif ordFound:
                dataVec.append(float(line.strip()))
            prevLine = line
        return np.array(dataVec, ndmin=2)
# ----------
    meshes = None
    sols   = None    

    for fp in fPaths:
        # Discarding MFEM index files
        if "mfem_root" not in fp and os.path.isdir(fp):
            meshesTmp = parseFileREMHOS(fp + "/mesh.000000")
            solsTmp   = parseFileREMHOS(fp + "/solution.000000")

            meshes = np.copy(meshesTmp) if meshes is None else np.concatenate((meshes, meshesTmp), axis=0)
            sols   = np.copy(solsTmp) if sols is None else np.concatenate((sols, solsTmp), axis=0)

    return [meshes, sols]

# ------------------------------------ CP2K Parser ------------------------------------

def loadDmpCP2K(fPaths):
    """
    A parser for CP2K data.
    
    This parser is able to handle Gaussian Cube files produced by
    CP2K. The specific metric represented in the given file is not
    relevant for the parser's operation.
  
    Parameters:
        fPaths (list[string]): Paths to files to be parsed
  
    Returns:
        list[np.array]: Matrices with parsed data
    """
    dmpMats = []

    for fp in fPaths:
        f = open(fp, 'r')
        numProps = 6
        numAtoms = 0
        lineIdx = 0
        tmpMat = []
        atomVec = []

        for line in f:
            splitLine = [el.strip() for el in line.split(" ") if el!=""]
            if lineIdx == 2:
                numAtoms = int(splitLine[0])
            elif lineIdx == 5:
                numProps = int(splitLine[0])
            elif lineIdx > (5 + numAtoms):
                atomVec = atomVec + parseFloats(splitLine)
                if len(atomVec) == numProps:
                    tmpMat.append(atomVec)
                    atomVec = []
            lineIdx = lineIdx + 1

        dmpMats.append(np.array(tmpMat, ndmin=2))

    return dmpMats

# ------------------------------------ QMCPACK Parser ------------------------------------

def loadDmpQMCPACK(fPaths):
    """
    A parser for QMCPack data.
    
    This function parses output data produced by QMCPack, specifically all
    system-level information within scalar.dat files. The parser discards
    columns that contain timing-related information (e.g., BlockCPU and 
    AcceptRatio) which is not part of the numerical solution. If multiple
    scalar.dat files are supplied, the respective columns are concatenated
    in sequence. The parser then returns one matrix (here 1D) for each 
    separate column in the files.
  
    Parameters:
        fPaths (list[string]): Paths to files to be parsed
  
    Returns:
        list[np.array]: Matrices with parsed data
    """
# ---------- Helper function
    def parseFileQMCPACK(fPathInt):
        f = open(fPathInt, 'r')
        tmpMat = []
        # Skipping first line
        f.readline()

        for line in f:
            splitLine = [el.strip() for el in line.split(" ") if el!=""]
            # Removing the block weight, block CPU and accept ratio properties
            tmpMat.append(parseFloats(splitLine[0:-3]))

        if len(tmpMat)==0:
            raise IOError("Empty QMCPack output file!")

        # In case block indexes are not sorted
        tmpMat.sort(key=lambda row: row[0])
        
        # Discarding index column
        return np.array(tmpMat, ndmin=2)[:,1:]
# ----------
    scalarDat = None

    for fp in fPaths:
        scalarTmp = parseFileQMCPACK(fp)
        scalarDat = np.copy(scalarTmp) if scalarDat is None else np.concatenate((scalarDat, scalarTmp), axis=0)

    # Separating single columns of the final matrix
    return [scalarDat[:, idx] for idx in range(scalarDat.shape[1])]

# ------------------------------------ LAMMPS Parser ------------------------------------

def loadDmpLAMMPS(fPaths):
    """
    A parser for LAMMPS data.
    
    This function parses basic LAMMPS dumps, containing coordinate 
    and velocity data for all atoms in a system over time. The parser
    returns one matrix for each dimension in coordinate and velocity
    vectors (so 6 in total). Within each matrix, the data for each 
    atom at each time step is stored.
  
    Parameters:
        fPaths (list[string]): Paths to files to be parsed
  
    Returns:
        list[np.array]: Matrices with parsed data
    """
    dmpMats = []

    for fp in fPaths:
        f = open(fp, 'r')
        numProps = 8
        numCoeffs = 6
        dmpDict = {}

        for line in f:
            splitLine = line.split(" ")
            # Hackish way to identify atom property lines
            if len(splitLine) == numProps:
                atomID = int(splitLine[0])
                atomVec = parseFloats(splitLine[2:])
                if not atomID in dmpDict:
                    dmpDict[atomID] = [atomVec]
                else:
                    dmpDict[atomID].append(atomVec)

        numAtoms = len(dmpDict.keys())
        numSteps = len(dmpDict[1])
        tmpMat = np.zeros((numAtoms, numSteps, numCoeffs))
        for k in range(numAtoms):
            for s in range(numSteps):
                for c in range(numCoeffs):
                    tmpMat[k][s][c] = dmpDict.get(k+1)[s][c]
        dmpMats = dmpMats + [np.array(tmpMat[:,:,idx], ndmin=2) for idx in range(numCoeffs)]

    return dmpMats
