#!/bin/bash

module load bowtie/1.3.1

MAPHELPER="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/software/hmmcopy_utils/build/util/mappability/generateMap.pl"
GENOME="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/crc-data/reference/hg38.fa"
WIG="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/crc-data/map/hg38.fa.map.wig"
BIGWIG="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/crc-data/map/hg38.fa.map.bw"
MAP="/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/software/hmmcopy_utils/build/bin/mapCounter"

$MAPHELPER -b $GENOME -o $BIGWIG

$MAPHELPER $GENOME -o $BIGWIG

$MAP -w 500000 $BIGWIG > "$WIG"

