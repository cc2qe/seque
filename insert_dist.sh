#!/bin/bash -e

if test -z "$1"
then
    echo "calculate the insert distribution from a BAM file"
    echo "usage: $0 <in.bam>"
    exit 1
fi

BAM=$1
BAMBASE=`basename $BAM .bam`
echo $BAMBASE
sambamba view -F "paired and mate_is_reverse_strand and not (unmapped or mate_is_unmapped or reverse_strand or secondary_alignment or duplicate)" $BAM \
    | awk '{ if ($7=="=" && $9>0 && $9<1000) print $9 }' \
    | sed -n '5000001,6000000p;6000000q' \
    | gzip -c \
    > $BAMBASE.insert.txt.gz

