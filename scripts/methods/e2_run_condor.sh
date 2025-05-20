#!/bin/bash

source /nfshomes/lmatalop/miniconda3/etc/profile.d/conda.sh
conda activate condor
source /etc/profile
module load gurobi

#PARAMETERS
NCELLS=$1
NCHARS=$2
NCLUSTERS=$3
k=$4
MISSING=$5
REPL=$6
NTHREADS=$7

#HYPERPARAMETERS
ESTIMATED_FP=0.0018
ESTIMATED_FN=0.16

MYMTHD="condor"
OUTPUT="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/simulation/n${NCELLS}_m${NCHARS}_p${NCLUSTERS}_k${k}_d${MISSING}/r${REPL}"
SIMDIR="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/simulation"
TOTAL="${SIMDIR}/n${NCELLS}_m${NCHARS}_p${NCLUSTERS}_k${k}_d${MISSING}/r${REPL}/gt_read_count.csv"
VARIANT="${SIMDIR}/n${NCELLS}_m${NCHARS}_p${NCLUSTERS}_k${k}_d${MISSING}/r${REPL}/gt_variant_count.csv"
MATRIX="${SIMDIR}/n${NCELLS}_m${NCHARS}_p${NCLUSTERS}_k${k}_d${MISSING}/r${REPL}/gt_character_matrix.csv"
CONDOR="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/software/ConDoR/src/condor.py"
PYTHON=$(which python)

echo "Running command:"
echo "$PYTHON $CONDOR -i $MATRIX.csv -a $ESTIMATED_FP -b $ESTIMATED_FN -k $k -r $TOTAL.csv -v $VARIANT.csv -o $OUTPUT"

MYTIME=$( { time $PYTHON $CONDOR \
    -i $MATRIX \
    -a $ESTIMATED_FP \
    -b $ESTIMATED_FN \
    -k $k \
    -r $TOTAL \
    -v $VARIANT \
    --threads $NTHREADS \
    -o ${OUTPUT}/${MYMTHD} &> ${OUTPUT}/${MYMTHD}.log ; } 2>&1 1>/dev/null)

uname -a > "${OUTPUT}/${MYMTHD}_node_info.csv"

MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
echo "$MYMTHD,$MYNODE,$MYSECS" > "${OUTPUT}/${MYMTHD}_runtime.csv"

