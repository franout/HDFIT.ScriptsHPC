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

from matplotlib.ticker import Locator
from enum import Enum
from bisect import bisect_right
import numpy as np
import logging as log
import re
import os

# Failure metric to be used in visualizations and other output
FI_METRIC = "AVF"

# Separators used to display instance trees on the shell
FI_SEPS = { "vert"  : "|",
            "int"   : "|-- ",
            "ext"   : "+-- ",
            "pad"   : "   ",
            "empty" : " ",
            "edge"  : "\"__ID__\" -> \"__ID2__\"\n",
            "shape" : "\"__ID__\" [label=\"__LABEL__\", shape=box]\n"}

# Enum used for indexing the different columns in an HDFIT CSV results file
class HDFIT(Enum):
    rank        = 0
    opCntExp    = 1
    opFi        = 2
    fiBit       = 3
    opCnt       = 4
    testPass    = 5
    testFail    = 6
    retCode     = 7
    time        = 8
    nrmse       = 9
    rtlError    = 10
    assignUUID  = 11
    moduleChain = 12
    moduleName  = 13

# ------------------------ Processing Netlist module info ------------------------

def parseModuleName(name, legacy=True):
    """
    Sanitizes a Netlist module name.
  
    This function removes all metadata within a Netlist module name string
    (e.g., the paramod, WIDTH and RESIDUE) fields, returning only the compact
    name of the module itself.
  
    Parameters:
        name (string): Full Netlist module name
        legacy (bool): If True, uses older VeriFI parsing conventions
  
    Returns:
        string: Sanitized version of name
    """
    nameL = [el.strip() for el in name.split("\\" if legacy else "\\\\") if el.strip()!=""]
    return nameL[0] if len(nameL)==1 else nameL[1]

def extractModules(fPath, legacy=True):
    """
    Parses module info from a VeriFI source file.
  
    This function is able to parse a VeriFI .cpp source file containing the
    static specification fo fault injection signals in a given Netlist (i.e.,
    detailing modules, signals contained therein and their features). It then
    returns a list of top-level module names that were identified, together 
    with the list of signal ID ranges associated with each. 
    
    BE CAREFUL: the values in sigs represent the first signal ID for each
    module, but these values are not guaranteed to be in sorted order.
  
    Parameters:
        fPath (string): Path to a .cpp VeriFI source file
        legacy (bool): If True, uses older VeriFI parsing conventions
  
    Returns:
        list[string]: Found module names
        list[int]:    Signal ID ranges for each module
    """
    f = open(fPath, "r")
    fText = f.read()
    f.close()

    mods = []
    sigs = []

    # Finding all module names
    for mod in re.finditer("{\n\s*\".+\",\n\s*{", fText):
        # Finding first signal within each module
        sig = re.search(",\n\s*[0-9]+,\n\s*},", fText[mod.span()[0]:-1])
        
        if sig is None:
            raise ValueError("Error encountered while extracting module information from %s!" % fPath)

        mods.append(parseModuleName(re.search("(?<=\").+(?=\")", mod.group(0)).group(0), legacy))
        sigs.append(int(re.search("[0-9]+(?=,)", sig.group(0)).group(0)))

    return np.asarray(mods), np.asarray(sigs)

def extractInstances(fPath, legacy=True):
    """
    Parses instance info from a VeriFI source file.
  
    This function parses a VeriFI .cpp source, similarly to extractModules(),
    and extract all information associated with module instances, which use their
    own UUID scheme. The function then returns a dictionary which maps each
    instance UUID (e.g., as found in the module instance chain CSV field) to
    the respective module name.
  
    Parameters:
        fPath (string): Path to a .cpp VeriFI source file
        legacy (bool): If True, uses older VeriFI parsing conventions
  
    Returns:
        dict{string : string}: Dictionary mapping UUIDs to instances
    """
    mods, sigs = extractModules(fPath, legacy)
    f = open(fPath, "r")
    fText = f.read()
    f.close()

    instDict = {}

    # Finding all instances
    for mod in re.finditer("(?<={)[0-9]+, [0-9]+(?=},)", fText):
        instTup = mod.group(0).split(", ")
        instDict[instTup[1]] = mods[int(instTup[0])]
    
    # Adding top element    
    instDict["1"] = "Top"
   
    return instDict

def searchModules(uuids, sigs):
    """
    Binary seach utility function.
  
    This function encapsulates a call to bisect_right in order to perform binary
    search on a given array. It is used to identify the modules to which the signals
    in a given input list belong.
  
    Parameters:
        uuids (list[int]): Signal IDs acting as keys
        sigs (list[int]): List of signal ID ranges over which to search
  
    Returns:
        list[int]: List of indexes of the module names to which the uuids belong
    """
    idxSigs = np.argsort(sigs)
    sortSigs = sigs[idxSigs]
    
    indexes = []
    for u in uuids:
        # Searching UUID within sorted list of module signal ranges
        uIdx = bisect_right(sortSigs, u) - 1
        # Mapping index to original unsorted list
        # A value of -1 implies no module was found
        indexes.append(idxSigs[uIdx] if uIdx >= 0 else -1)
       
    return indexes

def loadModuleData(uuids, fPath, legacy=True):
    """
    Utility function for loading module data.
  
    This utility function combines loading of module data from a VeriFI source
    (via extractModules) and mapping of a list of signal IDs to modules (via
    searchModules).
  
    Parameters:
        uuids (list[int]): List of signal IDs to be mapped to modules
        fPath (string): Path to a .cpp VeriFI source file
        legacy (bool): If True, uses older VeriFI parsing conventions
  
    Returns:
        list[string]: Module names associated with each uuids element
        list[string]: Found module names
        list[int]:    Signal ID ranges for each module
    """
    mods, sigs = extractModules(fPath, legacy)
    modsOut = [mods[el] if el>=0 else "NONE" for el in searchModules(uuids, sigs)]

    return np.asarray(modsOut), mods, sigs

def dumpModuleData(fPathIn, mods):
    """
    Augments an HDFIT CSV file with module name data
  
    This function will add a new column to an arbitrary HDFIT CSV results
    file, containing module name information for each injection signal, as
    computed by loadModuleData or in an equivalent way. The new CSV file can
    be distinguished from the original via the _modules suffix.
  
    Parameters:
        fPathIn (string): Path to an HDFIT CSV results file
        mods (list[string]): Module names associated with rows in fPathIn  
    """
    fIn  = open(fPathIn, "r")
    fOut = open(fPathIn[0:-4] + "_modules.csv", "w")
    
    header = fIn.readline()[0:-1] + ", module name\n"
    fOut.write(header)
    
    mIdx = 0
    for l in fIn:
        fOut.write(l[0:-1] + ", " + mods[mIdx] + "\n")
        mIdx = mIdx + 1
    
    fIn.close()
    fOut.close()

# ------------------------ Instance tree build functions ------------------------

def generateTree(instD, chains, failures, plot_rate, full_tree, tree_depth):
    """
    Generate a tree structure based on failure data
  
    This function generates a dictionary describing the RTL instance tree for
    a given hardware component, computing vulnerability metrics in the process.
    The function's output is a dictionary, in which every element is associated
    with a node in the tree: keys are represented by instance chains up to the 
    node's level (i.e., the path from the root to the node in the tree). Individual
    instances within the key may be identified by their raw UUID or by their 
    human-readable name according to the full_tree option.
    
    Each element's value is then represented by a tuple, in which the first field
    is a set containing the keys to the node's children, while the second and third
    field contain the number of failure occurrences and total injections respectively, 
    for the given node.
  
    Parameters:
        instD (dict{string : string}): Dictionary mapping UUIDs to instances,
            generated by the extractInstances function
        chains (list[string]): List of dash-separated UUID instance chains, one
            for each experiment run
        failures (list[bool]): List of failure occurrences associated with
            the supplied instance chains
        plot_rate (bool): If True, a failure rate value is generated for each node
            instead of the absolute number of failure occurrences
        full_tree (bool): If True, duplicated instances at the same tree level
            are not merged into one and the full tree is generated instead. Note
            that if this option is selected, the raw UUID instance chains are used 
            as keys to distinguish instances that have the same name
        tree_depth (int): Depth limit for the tree
    
    Returns:
        string: Key of the root node in the tree
        dict{string : (set(string), int, int)}): Dictionary describing
            the overall tree structure based on the instance chains
    """
    rootNode = instD["1"] if not full_tree else "1"
    instTree = {rootNode : [set(), 0, 0]}
    
    # Iterating over instance chains and updating nodes
    for i in range(len(chains)):
        # Checking for invalid instance chains
        if any(c=="0" for c in chains[i].split("-")):
            log.warning("Invalid instance chain %s detected. Skipping..." % chains[i])
            continue
        iChain = [instD[c] if not full_tree else c for c in chains[i].split("-")]
        iDepth = int(np.min((tree_depth+1, len(iChain))))
        node   = iChain[0]
        for j in range(1, iDepth):
            node   = "-".join(iChain[0:j+1])
            parent = "-".join(iChain[0:j])
            instTree = updateNode(instTree, parent, node, failures[i])
        instTree = updateNode(instTree, node, None, failures[i])
    
    # Adjusting failure occurrences into failure rate values if required
    if plot_rate:
        for k in instTree.keys():
            instTree[k][1] = instTree[k][1]/instTree[k][2]*100 if instTree[k][2]!=0 else 0
            
    return rootNode, instTree

def updateNode(treeDict, key, child, fail_occ):
    """
    Internal function for updating instance tree elements
  
    This function is meant for internal use only. It updates the tree dictionary
    maintained by the generateTree function, adding new information about nodes.
  
    Parameters:
        treeDict (dict{string : (set(string), int, int)}): Dictionary describing
            the tree structure, created by the generateTree function
        key (string): Key of node to be updated
        child (string): Key of child to be added to node
        fail_occ (bool): If True, adds a new failure occurrence to the node
    
    Returns:
        dict{string : (set(string), int, int)}: Updated tree dict structure
    """
    if key not in treeDict:
        treeDict[key] = [set(), 0, 0]
    if child is not None:
        treeDict[key][0].add(child)
    treeDict[key][1] = treeDict[key][1] + fail_occ
    treeDict[key][2] = treeDict[key][2] + 1
    return treeDict

def displayTree(key, instTree, instDict, plot_rate, lc=True, padPre=""):
    """
    Display an HDFIT instance tree on the shell
  
    This function is able to display an HDFIT instance tree built using
    the generateTree function on the shell, optionally creating a Dot
    code representation that can be visualized through any external
    Graphviz engine. This function is recursive: the lc and padPre 
    arguments are meant to steer the internal recursion and are not
    to be set externally.
  
    Parameters:
        key (string): Key of the root node from which to start
        instTree (dict{string : (set(string), int, int)}): Dictionary describing
            the tree structure, created by the generateTree function
        instDict (dict{string : string}): Dictionary mapping UUIDs to instances,
            generated by the extractInstances function
        plot_rate (bool): If True, the failure rate is shown, else the occurrences
    
    Returns:
        string: String with representation for Dot visualization
    """
    nodeName  = formatTreeNode(key, instTree, instDict, plot_rate)
    # Determining padding string for this node's children
    padPost = padPre + (FI_SEPS["vert"] if not lc else FI_SEPS["empty"]) + FI_SEPS["pad"]
    # Adding new shape to graphviz dot string
    shapes  = FI_SEPS["shape"].replace("__ID__", key).replace("__LABEL__", nodeName)
    edges   = ""
    
    # Display string for this node
    print(padPre + (FI_SEPS["int"] if not lc else FI_SEPS["ext"]) + nodeName)
    
    # Constructing ordered list of child - failure rate pairs
    childList = [[el, instTree[el][1]] for el in instTree[key][0]]
    childList.sort(key = lambda x : x[1], reverse = True)
    # Iterating over children
    for idx in range(len(childList)): 
        c   = childList[idx]
        # Iterating over children
        s, e  = displayTree(c[0], instTree, instDict, plot_rate, idx >= len(childList) - 1, padPost)
        # Adding new shapes and edges
        shapes  = shapes + s
        edges   = edges + FI_SEPS["edge"].replace("__ID__", key).replace("__ID2__", c[0]) + e
    return formatGraphviz(shapes, edges) if padPre=="" else (shapes, edges)
    
def formatTreeNode(key, instTree, instDict, plot_rate):
    """
    Generates a string for showing a single instance tree node
  
    This utility function generates a string for a given node in the instance
    tree, showing its name coupled with the associated failure rate or absolute 
    number of occurrences.
  
    Parameters:
        key (string): Key of the target node
        instTree (dict{string : (set(string), int, int)}): Dictionary describing
            the tree structure, created by the generateTree function
        instDict (dict{string : string}): Dictionary mapping UUIDs to instances,
            generated by the extractInstances function
        plot_rate (bool): If True, the failure rate is shown, else the occurrences
    
    Returns:
        string: String with information for the given node
    """
    frStr    = str(instTree[key][1]) if not plot_rate else '{:.4f}'.format(instTree[key][1]) + "%"
    nodeName = (key.split("-")[-1] if instDict is None else instDict[key.split("-")[-1]]) + ": " + frStr
    return nodeName
    
def formatGraphviz(shapes, edges):
    """
    Generates a Dot representation for a given instance tree
  
    This function takes as input the list of shapes and edges created
    within the displayTree function, and outputs a Dot representation
    that can be displayed via any Graphviz engine.
  
    Parameters:
        shapes (string): String containing Dot shapes
        edges (string): String containing Dot edges
    
    Returns:
        string: String with representation for Dot visualization
    """
    dotStr = "digraph tree {\n"
    dotStr = dotStr + "rankdir=\"LR\"\n" 
    dotStr = dotStr + shapes + "\n" + edges + "}\n"
    return dotStr

# ------------------------ Data loading and processing ------------------------

def loadData(fault_log):
    """
    Reads fault data from a CSV file.
  
    This function parses an HDFIT CSV results file and returns a dictionary,
    where each element is a list corresponding to a column in the file. 
    Elements in the dictionary are indexed via the HDFIT enum class. Note
    that not all of the column available in CSV files are currently parsed,
    but only those that are actively needed by other plotting functions.
  
    Parameters:
        fault_log (string): Path to an HDFIT CSV results file
  
    Returns:
        dict: Dictionary with parsed HDFIT data
    """
    fDict = {}
    # Read failure results from the log file
    fDict[HDFIT.opCntExp] = np.genfromtxt(fault_log, skip_header=1, delimiter=',', usecols=1, dtype=np.int64)
    fDict[HDFIT.opFi] = np.genfromtxt(fault_log, skip_header=1, delimiter=',', usecols=2, dtype=np.int64)
    fDict[HDFIT.fiBit] = np.genfromtxt(fault_log, skip_header=1, delimiter=',', usecols=3, dtype=int)
    fDict[HDFIT.opCnt] = np.genfromtxt(fault_log, skip_header=1, delimiter=',', usecols=4, dtype=np.int64)
    fDict[HDFIT.testFail] = np.genfromtxt(fault_log, skip_header=1, delimiter=',', usecols=6, dtype=int)
    # Read SDE errors (if available)
    try:
        fDict[HDFIT.nrmse] = np.genfromtxt(fault_log, skip_header=1, delimiter=',', usecols=9, dtype=float)
    except ValueError:
        log.warning("File %s does not appear to contain NRMSE error data." % fault_log)
        fDict[HDFIT.nrmse] = np.zeros(fDict[HDFIT.testFail].shape)
    # Read RTL error status (if available)
    try:
        fDict[HDFIT.rtlError] = np.genfromtxt(fault_log, skip_header=1, delimiter=',', usecols=10, dtype=int)
        fDict[HDFIT.assignUUID] = np.genfromtxt(fault_log, skip_header=1, delimiter=',', usecols=11, dtype=int)
        fDict[HDFIT.moduleChain] = np.genfromtxt(fault_log, skip_header=1, delimiter=',', usecols=12, dtype=str)
    except ValueError:
        log.warning("File %s does not appear to contain RTL-related data." % fault_log)
        fDict[HDFIT.rtlError] = np.zeros(fDict[HDFIT.testFail].shape)
        fDict[HDFIT.assignUUID] = np.zeros(fDict[HDFIT.testFail].shape)
        fDict[HDFIT.moduleChain] = np.zeros(fDict[HDFIT.testFail].shape)

    return fDict

def filterData(dataDict, bit_pos, exclude_uoe, exclude_rtl):
    """
    Produces a mask for excluding runs from an HDFIT experiment.
  
    This function computes a "mask" list starting from a data
    dictionary loaded via loadData, associating True to any 
    given element if the corresponding experiment run should be
    excluded from failure rate computations, and False otherwise. This
    is determined based on a series of user-configurable settings.
  
    Parameters:
        dataDict (dict): Dictionary containing HDFIT experiment data
        bit_pos (int): Exclude runs with faults injected in bits higher than bit_pos
        exclude_uoe (bool): If True, exclude all runs resulting in observable failures
        exclude_rtl (bool): If True, exclude all runs where RTL errors were triggered
  
    Returns:
        list[bool]: Mask list dictating which experiment runs are to be ignored
    """
    failure   = dataDict[HDFIT.testFail]
    fault     = dataDict[HDFIT.fiBit]
    rtlerrors = dataDict[HDFIT.rtlError]
    max_row = len(failure)
    
    if max_row==0:
        raise ValueError("The input HDFIT dataset is empty!")
    
    # Computing runs to be excluded from the failure rate computation
    ex_runs = [fault[i]>=bit_pos or (failure[i] and exclude_uoe) or (rtlerrors[i] and exclude_rtl) for i in range(max_row)] 
    return ex_runs

def getDataPartialProtection(dataDict, bit_pos, bins, exclude_uoe, exclude_rtl):
    """
    Produce an NRMSE vs failure rate (i.e., PVF or AVF) curve.
  
    This function produces values for a NRMSE vs failure rate curve according to an
    arbitrary set of bin values, as well as flags to determine which runs
    should be excluded using the filterData function. If a given run is 
    to be excluded, it will be considered successful in the computation.
  
    Parameters:
        dataDict (dict): Dictionary containing HDFIT experiment data
        bit_pos (int): Exclude runs with faults injected in bits higher than bit_pos
        bins (list[float]): List of bin edges for which failure rate values are computed
        exclude_uoe (bool): If True, exclude all runs resulting in observable failures
        exclude_rtl (bool): If True, exclude all runs where RTL errors were triggered
  
    Returns:
        list[float]: List of failure rate values for the points specified by bins
    """
    failure   = dataDict[HDFIT.testFail]
    fault     = dataDict[HDFIT.fiBit]
    errors    = dataDict[HDFIT.nrmse]
    rtlErrors = dataDict[HDFIT.rtlError]
    max_row = len(failure)
    
    if max_row==0:
        raise ValueError("The input HDFIT dataset is empty!")
    
    # List storing the final curve
    vals = np.zeros(bins.shape)
    # Filtering out runs
    ex_runs = filterData(dataDict, bit_pos, exclude_uoe, exclude_rtl)
    # Computing list of effective failures (errors and SDEs) after filtering out runs
    temp_errors = [(errors[i]>0 or failure[i]>0) and not ex_runs[i] for i in range(max_row)]
    
    for idx in range(len(bins)):
        mask = (errors > bins[idx])
        vals[idx] = np.sum(temp_errors, where=mask) / max_row
    
    return vals

# --------------------  PVF or AVF statistics computation  --------------------

def printFailureStatistics(dataDict, bit_pos, bit_depth):
    """
    Print a variety of failure rate (i.e., PVF or AVF) statistics.
  
    This function takes as input a dictionary containing HDFIT experiment data,
    and will compute a breakdown of the main failure metrics (e.g., UOEs and SDEs,
    faults caught by RTL protection, and so on). A "partial protection" setting
    can also be applied, ignoring all failures occurring in bits over a certain order.
  
    Parameters:
        dataDict (dict): Dictionary containing HDFIT experiment data
        bit_pos (int): Exclude runs with faults injected in bits higher than bit_pos
        bit_depth (int): Bit precision currently being used
    """    
    failure   = dataDict[HDFIT.testFail]
    fault     = dataDict[HDFIT.fiBit]
    errors    = dataDict[HDFIT.nrmse]
    rtlerrors = dataDict[HDFIT.rtlError]
    max_row   = len(failure)
    
    if max_row==0:
        raise ValueError("The input HDFIT dataset is empty!")
    
    uoeFR = np.sum([failure[i] and fault[i]<bit_pos for i in range(max_row)]) * 100 / max_row
    sdeFR = np.sum([errors[i] > 0 and not failure[i] and fault[i]<bit_pos for i in range(max_row)]) * 100 / max_row
    
    # Computing number of runs where RTL protection was triggered
    rtlRuns = [rtlerrors[i] and fault[i]<bit_pos for i in range(max_row)]
    rtlNum  = np.sum(rtlRuns) * 100 / max_row
    rtlFR_uoe = np.sum([rtlRuns[i] and failure[i] for i in range(max_row)]) * 100 / max_row
    rtlFR_sde = np.sum([rtlRuns[i] and errors[i] > 0 and not failure[i] for i in range(max_row)]) * 100 / max_row
    rtlFR_fp  = rtlNum - rtlFR_uoe - rtlFR_sde
    # Adjusting overall statistics according to RTL protection
    uoeFR      = uoeFR - rtlFR_uoe
    sdeFR      = sdeFR - rtlFR_sde
    overallFR  = uoeFR + sdeFR
    successRuns = 100 - overallFR - rtlNum
    
    if bit_pos < bit_depth:
        log.info("-- %s with protection up to bit %d: %f%%" % (FI_METRIC, bit_pos, overallFR))
    else:
        log.info("-- %s with no protection: %f%%" % (FI_METRIC, overallFR))
    log.info("---- UOE %s    : %f%%" % (FI_METRIC, uoeFR))
    log.info("---- SDE %s    : %f%%" % (FI_METRIC, sdeFR))
    log.info("---- Success    : %f%%" % successRuns)
    if rtlNum > 0:
        log.info("---- RTL Errors : %f%%" % rtlNum)
        log.info("-------- Leading to UOE     : %f%%" % rtlFR_uoe)
        log.info("-------- Leading to SDE     : %f%%" % rtlFR_sde)
        log.info("-------- Leading to Success : %f%%" % rtlFR_fp)

# ------------------------     Miscellaneous code      ------------------------
 
def seekRightBound(errors):
    """
    Returns the 95th percentile of the input list.
  
    This function is used to determine a suitable right bound when 
    plotting an error distribution, assuming the left bound is 0. 
    The 95th percentile usually works well in practice.
  
    Parameters:
        errors (list[float]): List of error values
  
    Returns:
        float: 95th percentile of errors list
    """
    errors = [ val for val in errors if val>0 and val!=np.inf ]
    # Picking the largest statistically relevant error as right bound
    rb = np.percentile(errors, 95)
    return rb
    
def getFileBaseName(name):
    """
    Removes date and extension from an HDFIT CSV file name.
    
    Parameters:
        name (string): Name of an HDFIT CSV file
  
    Returns:
        string: Sanitized version of name
    """
    baseName = os.path.basename(name)
    match = re.search("[-_.][0-9]+-[0-9]+-[0-9]+.csv", baseName)
    return baseName.replace(match.group(0), "") if match is not None else baseName[0:-4]

def percent(arr):
    """
    Calculates an array's mean as a percentage.
    
    This utility function is used as an estimator in
    some seaborn barplot calls, to plot failure rate
    bars (which are expressed as percentages). No
    normalization is performed, and the range of data
    is assumed to be [0, 1].
    
    Parameters:
        arr (np.array): Data to be processed
  
    Returns:
        float: Percentage mean of the array
    """    
    arr = arr.astype(float)    
    return np.mean(arr) * 100

class MinorSymLogLocator(Locator):
    """
    This is a utility class to produce ticks in symlog scale.
    """
    def __init__(self, linthresh):
        self.linthresh = linthresh

    def __call__(self):
        majlocs = self.axis.get_majorticklocs()
        minlocs = []
        for i in range(1, len(majlocs)):
            majstep = majlocs[i] - majlocs[i-1]
            ndivs = 10 if abs(majlocs[i-1] + majstep/2) < self.linthresh else 9
            minstep = majstep / ndivs
            locs = np.arange(majlocs[i-1], majlocs[i], minstep)[1:]
            minlocs.extend(locs)

        return self.raise_if_exceeds(np.array(minlocs))

    def tick_values(self, vmin, vmax):
        raise NotImplementedError("")

def HDFIT_intList(s):
    """
    Splits a comma-separated string into a list of integers.

    Parameters:
        s (string): String to be separated

    Returns:
        list[int]: List of integers extracted from s
    """
    ls = [int(el.strip()) for el in s.split(",") if el.strip()!=""]
    if any(el < 0 for el in ls):
        raise ValueError("List %s contains one or more invalid values. Plase check your input." % s)
    return ls

def HDFIT_intPos(s):
    """
    Converts a string into a (positive) integer.

    Parameters:
        s (string): String to be converted

    Returns:
        int: Integer representation of s
    """
    i = int(s)
    if i <= 0:
        raise ValueError("Value %s is invalid. Please check your input." % s)
    return i
