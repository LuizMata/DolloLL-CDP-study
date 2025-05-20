#!/bin/bash

source /nfshomes/lmatalop/miniconda3/etc/profile.d/conda.sh

conda activate condor

echo "Current Conda environment: $CONDA_DEFAULT_ENV"

#PARAMETERS
NCELLS=$1
NCHARS=$2
NCLUSTERS=$3
k=$4
MISSING=$5
REPL=$6

#HYPERPARAMETERS
FP=0.001
FN=0.001
LAMBDA=0.8
VAF_THRESHOLD=0.1
READ_THRESHOLD=5
ESTIMATED_FP=0.0018
ESTIMATED_FN=0.16

#DIRECTORIES

SIMDIR="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/simulation"
OUTPREFIX="${SIMDIR}/n${NCELLS}_m${NCHARS}_p${NCLUSTERS}_k${k}_d${MISSING}/r${REPL}"
SIMSCRIPT="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/util/simulation_reads_simple.py"
PYTHON=$(which python)

$PYTHON $SIMSCRIPT \
-s $REPL \
-o $OUTPREFIX \
-n $NCELLS \
-m $NCHARS \
-d $MISSING \
-p $NCLUSTERS \
-a $FP \
-b $FN \
-k $k \
-l $LAMBDA -v \
&> "${OUTPREFIX}/simulation.log" 2>&1 1> "${OUTPREFIX}/simulation.err"

