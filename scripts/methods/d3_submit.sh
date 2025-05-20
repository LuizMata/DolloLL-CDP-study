#!/bin/bash

#PARAMETERS
NCELLS=( 25 50 100 )    
NCLUSTERS=( 5 ) 
MISSING=( 0 )
K=( 3 )
REPLS=( $(seq 1 2) )

for NCELL in ${NCELLS[@]}; do
    for NCLUSTER in ${NCLUSTERS[@]}; do
        for k in ${K[@]}; do
            for d in ${MISSING[@]}; do
                for REPL in ${REPLS[@]}; do

echo "executing Dollo-CDPLLRC - cells:$NCELL, clusters:$NCLUSTER, k:$k, missing:$d, rep:$REPL..."

sbatch \
    --job-name="Dollo-CDPLLRC.$NCELL.$NCLUSTER.$k.$d.$REPL" \
    --export=NCELL="$NCELL",NCLUSTER="$NCLUSTER",k="$k",MISSING="$d",REPL="$REPL" \
d1_drive.sbatch

	        done
            done
        done
    done
done


