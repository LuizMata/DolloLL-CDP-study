#!/bin/bash

source /nfshomes/lmatalop/miniconda3/etc/profile.d/conda.sh
conda activate condor
source /etc/profile
module load gurobi

NCELLS=$1
NCHARS=$2
NCLUSTERS=$3
k=$4
MISSING=$5
REPL=$6
METHOD=$7

OUTPUT="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/simulation/n${NCELLS}_m${NCHARS}_p${NCLUSTERS}_k${k}_d${MISSING}/r${REPL}/${METHOD}_evaluation_contracted"
SIMDIR="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/simulation"
EVALUATE="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/util/evaluate_reads_simple.py"
PYTHON=$(which python)

$PYTHON $EVALUATE \
	-s $REPL \
	--method $METHOD \
	-o $OUTPUT.csv \
	-n $NCELLS \
	-m $NCHARS \
	-d $MISSING \
	-p $NCLUSTERS \
	-k $k \
	--dir $SIMDIR \
	> $OUTPUT.out 2> $OUTPUT.err
