#!/bin/bash

source /nfshomes/lmatalop/miniconda3/etc/profile.d/conda.sh

conda activate condor

NCELLS=$1
NCHARS=$2
NCLUSTERS=$3
k=$4
MISSING=$5
REPL=$6

MYMTHD="dollocdpLLRC"

OUTPUT="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/simulation/n${NCELLS}_m${NCHARS}_p${NCLUSTERS}_k${k}_d${MISSING}/r${REPL}"
SIMDIR="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/simulation"
TOTAL="${SIMDIR}/n${NCELLS}_m${NCHARS}_p${NCLUSTERS}_k${k}_d${MISSING}/r${REPL}/read_count.csv"
VARIANT="${SIMDIR}/n${NCELLS}_m${NCHARS}_p${NCLUSTERS}_k${k}_d${MISSING}/r${REPL}/variant_count.csv"
BED="${SIMDIR}/n${NCELLS}_m${NCHARS}_p${NCLUSTERS}_k${k}_d${MISSING}/r${REPL}/locations.bed"
DOLLOCDPLL="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/software/Dollo-CDPLLRC/src/dollo-cdp"

MYTIME=$(time ($DOLLOCDPLLRC \
    -t1 $TOTAL \
    -v $VARIANT \
    -l $BED \
    -d 3 \
    -g O \
    -o ${OUTPUT}/${MYMTHD}.tre \
    &> ${OUTPUT}/${MYMTHD}.log) 2>&1 1>/dev/null)

uname -a > "${OUTPUT}/${MYMTHD}_node_info.csv"

MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
echo "$MYMTHD,$MYNODE,$MYSECS" > "${OUTPUT}/${MYMTHD}_runtime.csv"

