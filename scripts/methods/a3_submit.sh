#!/bin/bash

#PARAMETERS
NCELLS=( 25 50 100 )    
NCLUSTERS=( 3 5 ) 
MISSING=( 0 0.1 )
K=( 1 2 3 )
REPLS=( $(seq 1 5) )

for NCELL in ${NCELLS[@]}; do
    for NCLUSTER in ${NCLUSTERS[@]}; do
        for k in ${K[@]}; do
            for d in ${MISSING[@]}; do
                for REPL in ${REPLS[@]}; do

echo "simulating - cells:$NCELL, clusters:$NCLUSTER, k:$k, missing:$d, rep:$REPL..."

sbatch \
    --job-name="sim-cond-a1.$NCELL.$NCLUSTER.$k.$d.$REPL" \
    --export=NCELL="$NCELL",NCLUSTER="$NCLUSTER",k="$k",MISSING="$d",REPL="$REPL" \
a1_drive.sbatch

	        done
            done
        done
    done
done

