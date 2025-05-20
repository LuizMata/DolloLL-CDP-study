#!/bin/bash

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

#MAIN="n${NCELL}_m${NCELL}_p${NCLUSTER}_k${k}_d${d}"
#mkdir $MAIN

SUB="n${NCELL}_m${NCELL}_p${NCLUSTER}_k${k}_d${d}/r${REPL}"

mkdir $SUB

	        done
            done
        done
    done
done

