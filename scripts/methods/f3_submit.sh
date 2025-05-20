#!/bin/bash

#PARAMETERS
NCELLS=( 25 50 100 )    
NCLUSTERS=( 3 5 ) 
MISSING=( 0 0.1 )
K=( 1 2 3 )
REPLS=( $(seq 1 5) )
METHODS=( "condor" "dollocdp" )

for METHOD in ${METHODS[@]}; do
    for NCELL in ${NCELLS[@]}; do
        for NCLUSTER in ${NCLUSTERS[@]}; do
            for k in ${K[@]}; do
                for d in ${MISSING[@]}; do
                    for REPL in ${REPLS[@]}; do

echo "Running Evaluation..."

sbatch \
    --job-name="Evaluation.$NCELL.$NCLUSTER.$k.$d.$REPL" \
    --export=NCELL="$NCELL",NCLUSTER="$NCLUSTER",k="$k",MISSING="$d",REPL="$REPL",METHOD="$METHOD" \
f1_drive.sbatch

		    done
	        done
            done
        done
    done
done


