#!/bin/bash

#PARAMETERS
NCELLS=$1
NCHARS=$2
NCLUSTERS=$3
k=$4
MISSING=$5
REPL=$6
METHOD=$7

#DIRECTORIES

SIMDIR="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/simulation"
INTREE="${SIMDIR}/n${NCELLS}_m${NCHARS}_p${NCLUSTERS}_k${k}_d${MISSING}/r${REPL}/${METHOD}_tree.newick"
INCHAR="${SIMDIR}/n${NCELLS}_m${NCHARS}_p${NCLUSTERS}_k${k}_d${MISSING}/r${REPL}/gt_nexus2.nex"
CONTRACTION="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/util/tree_util/edgeContractor"

$CONTRACTION $INTREE $INCHAR

