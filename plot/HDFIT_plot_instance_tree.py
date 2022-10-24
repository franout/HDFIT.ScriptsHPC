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

from HDFIT_plot_utils import *
import sys, argparse
import logging as log
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def main():
    log.basicConfig(format="[ %(levelname)s ] %(message)s", level=log.INFO, stream=sys.stdout)

    # Input parsing
    parser = argparse.ArgumentParser(description="HDFIT HW module instance tree plotting script")
    parser.add_argument("path", action="store", type=str, help="Path to HDFIT CSV file")
    parser.add_argument("signals", action="store", type=str, help="C++ file with FI signal spec")
    parser.add_argument("-b", action="store", dest="bit_depth", type=HDFIT_intPos, default=sys.maxsize, help="Maximum fault injection bit")
    parser.add_argument("-d", action="store", dest="tree_depth", type=HDFIT_intPos, default=sys.maxsize, help="Maximum depth of the tree")
    parser.add_argument("--uoe", action="store_true", dest="exclude_uoe", help="Exclude runs resulting in UOE")
    parser.add_argument("--rtl", action="store_false", dest="exclude_rtl", help="Include runs containing RTL errors")
    parser.add_argument("--rate", action="store_true", dest="plot_rate", help="Display failure rate instead of occurrences")
    parser.add_argument("--full", action="store_true", dest="full_tree", help="Display full instance tree")
    parser.add_argument("--viz", action="store_true", dest="dmp_graphviz", help="Dump Dot representation for Graphviz rendering")
    args = parser.parse_args()
    
    # Getting parsed arguments
    bit_depth    = args.bit_depth
    tree_depth   = args.tree_depth
    sig_file     = args.signals
    exclude_uoe  = args.exclude_uoe
    exclude_rtl  = args.exclude_rtl
    plot_rate    = args.plot_rate
    full_tree    = args.full_tree
    dmp_graphviz = args.dmp_graphviz
    fault_log    = args.path

    # Loading and unpacking data
    dataDict  = loadData(fault_log)
    uuid      = dataDict[HDFIT.assignUUID]
    failure   = dataDict[HDFIT.testFail]
    errors    = dataDict[HDFIT.nrmse]
    chains    = dataDict[HDFIT.moduleChain]
    
    max_row = len(failure)
    log.info(f"{max_row} faults loaded from the fault log file")
    
    # Filtering out runs
    ex_runs = filterData(dataDict, bit_depth, exclude_uoe, exclude_rtl)
    # Combining hard program failures with SDEs
    agg_failures = np.asarray([(failure[idx]>0 or errors[idx]>0) and not ex_runs[idx] for idx in range(max_row)])
    
    # Constructing tree and failure rate values
    instD = extractInstances(sig_file)
    rootNode, instTree = generateTree(instD, chains, agg_failures, plot_rate, full_tree, tree_depth)
    
    # Displaying final tree structure
    log.info("Module instance tree:")
    gviz = displayTree(rootNode, instTree, None if not full_tree else instD, plot_rate)
    
    # Dumping Graphviz dot representation
    if dmp_graphviz:
        f = open(fault_log[0:-4] + "_instance_tree.dot", "w")
        f.write(gviz)
        f.close()
 
if __name__ == '__main__':
    sys.exit(main() or 0)
