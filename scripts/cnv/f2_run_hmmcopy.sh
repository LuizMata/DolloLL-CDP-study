#!/bin/bash

source /nfshomes/lmatalop/miniconda3/etc/profile.d/conda.sh
conda activate visuals

CELL=$1

Rscript hmmcopy.R $CELL

