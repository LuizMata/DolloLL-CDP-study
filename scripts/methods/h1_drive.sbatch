#!/bin/bash

# https://slurm.schedmd.com/srun.html
# https://rc-docs.northeastern.edu/en/latest/using-discovery/srun.html

#SBATCH --qos=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64gb
#SBATCH --mem-bind=local
#SBATCH --account=cbcb
#SBATCH --partition=cbcb
#SBATCH --time=00:30:00

srun \
    --nodes 1 \
    --ntasks 1 \
    --cpus-per-task 1 \
    --mem=64gb \
    ./h2_edge_contraction.sh $NCELL $NCELL $NCLUSTER $k $MISSING $REPL $METHOD
wait

