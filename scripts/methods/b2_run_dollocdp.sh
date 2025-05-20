#!/bin/bash

#PARAMETERS
NCELLS=$1
NCHARS=$2
NCLUSTERS=$3
k=$4
MISSING=$5
REPL=$6

MYMTHD="dollocdp"

SIMDIR="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/simulation"
OUTPUT="${SIMDIR}/n${NCELLS}_m${NCHARS}_p${NCLUSTERS}_k${k}_d${MISSING}/r${REPL}"
MATRIX="${SIMDIR}/n${NCELLS}_m${NCHARS}_p${NCLUSTERS}_k${k}_d${MISSING}/r${REPL}/gt_nexus.nex"
DOLLOCDP="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/software/Dollo-CDP/src/dollo-cdp"

MYTIME=$(time ($DOLLOCDP \
    -i $MATRIX \
    -g O \
    -o ${OUTPUT}/${MYMTHD}.tre \
    &> ${OUTPUT}/${MYMTHD}.log) 2>&1 1>/dev/null)

uname -a > "${OUTPUT}/${MYMTHD}_node_info.csv"

MYNODE=$( uname -a | sed 's/\./ /g' | awk '{print $2}' )

MYSECS="$(echo $MYTIME | awk '{print $2","$4","$6}')"
echo "$MYMTHD,$MYNODE,$MYSECS" > "${OUTPUT}/${MYMTHD}_runtime.csv"



