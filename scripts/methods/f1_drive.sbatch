#!/bin/bash

# https://slurm.schedmd.com/srun.html
# https://rc-docs.northeastern.edu/en/latest/using-discovery/srun.html

#SBATCH --qos=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --mem-bind=local
#SBATCH --account=cbcb
#SBATCH --partition=cbcb
#SBATCH --time=1:00:00
#SBATCH --constraint=EPYC-7313

srun \
    --nodes 1 \
    --ntasks 1 \
    --cpus-per-task 1 \
    --mem=8gb \
    ./f2_run_evaluation.sh $NCELL $NCELL $NCLUSTER $k $MISSING $REPL $METHOD

wait

