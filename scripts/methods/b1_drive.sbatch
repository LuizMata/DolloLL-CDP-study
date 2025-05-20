#!/bin/bash

# https://slurm.schedmd.com/srun.html
# https://rc-docs.northeastern.edu/en/latest/using-discovery/srun.html

#SBATCH --qos=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1024gb
#SBATCH --mem-bind=local
#SBATCH --exclusive
#SBATCH --account=cbcb
#SBATCH --partition=cbcb
#SBATCH --time=18:00:00
#SBATCH --constraint=EPYC-7313

srun \
    --nodes 1 \
    --ntasks 1 \
    --cpus-per-task 1 \
    --mem=1024gb \
    ./b2_run_dollocdp.sh $NCELL $NCELL $NCLUSTER $k $MISSING $REPL
wait
