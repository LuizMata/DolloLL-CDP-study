#!/bin/bash

CELLS=()
while IFS= read -r line; do
    CELLS+=("$line")
done < ../../crc-data/cell-fastq/cells.txt

for CELL in ${CELLS[@]}; do

sbatch \
	--job-name="hmmcopy" \
	--export=CELL="$CELL" \
	f1_drive.sbatch

done
	
