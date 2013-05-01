#!/bin/bash

BASE_CHRS="\
chr1 \
chr2 \
chr3 \
chr4 \
chr5 \
chr6 \
chr7 \
chr8 \
chr9 \
chr10 \
chr11 \
chr12 \
chr13 \
chr14 \
chr15 \
chr16 \
chr17 \
chr18 \
chr19 \
chr20 \
chr21 \
chr22 \
chr23 \
chr24 \
chr25 \
chr26 \
chr27 \
chr28 \
chr29 \
chrX \
"

TRANSGENE="HD_transgene"
SUM=0
NUM_BASES=0

while IFS=$'\t' read -r -a myBed
do
    chrom=${myBed[0]}
    start=${myBed[1]}
    end=${myBed[2]}
    depth=${myBed[3]}

    if [[ "$BASE_CHRS" == *"$chrom"* ]]
    then
	let "NUM_BASES += $(( $end - $start ))"
	let "SUM += $(( $depth * ($end - $start) ))"
    fi
    
done < $1

echo $SUM, $NUM_BASES