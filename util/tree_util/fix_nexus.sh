#!/bin/bash

#PARAMETERS
NCELLS=( 25 50 100 )    
NCLUSTERS=( 3 5 ) 
MISSING=( 0 0.1 )
K=( 1 2 3 )
REPLS=( $(seq 1 5) )
SIMDIR="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/simulation"

for NCELL in ${NCELLS[@]}; do
    for NCLUSTER in ${NCLUSTERS[@]}; do
        for k in ${K[@]}; do
            for d in ${MISSING[@]}; do
                for REPL in ${REPLS[@]}; do

INPUT="${SIMDIR}/n${NCELL}_m${NCELL}_p${NCLUSTER}_k${k}_d${d}/r${REPL}/gt_nexus.nex"
OUTPUT="${SIMDIR}/n${NCELL}_m${NCELL}_p${NCLUSTER}_k${k}_d${d}/r${REPL}/gt_nexus2.nex"

sed -e "s/ntax=$((NCELL+1))/ntax=${NCELL}/g" -e '/^O/d' "${INPUT}" > "${OUTPUT}"

	        done
            done
        done
    done
done


