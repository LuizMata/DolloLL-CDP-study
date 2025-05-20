#!/bin/bash

# https://slurm.schedmd.com/srun.html
# https://rc-docs.northeastern.edu/en/latest/using-discovery/srun.html

#SBATCH --qos=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256gb
#SBATCH --mem-bind=local
#SBATCH --account=cbcb
#SBATCH --partition=cbcb
#SBATCH --time=12:00:00
#SBATCH --constraint=EPYC-7313


srun \
    --nodes 1 \
    --ntasks 1 \
    --cpus-per-task 1 \
    --mem=256gb \
    ./f2_run_hmmcopy.sh $CELL 

wait

