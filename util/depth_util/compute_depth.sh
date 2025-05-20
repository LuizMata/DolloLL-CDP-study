#!/bin/bash

source /etc/profile
module load samtools

INPUT="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/crc-data/alignments/SRR3472835-p8.sorted.bam"
OUTPUT="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/util/average_depth.txt"

samtools depth -a $INPUT | awk '{sum+=$3} END {if (NR>0) print sum / NR; else print "0"}' > "$OUTPUT"

