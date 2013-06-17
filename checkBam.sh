#!/bin/bash

############################################################
#  Program: checkBam
#  Author : Colby Chiang
#  Description: Checks a bam file for truncation
############################################################


## BEGIN SCRIPT
usage()
{
    cat << EOF

usage: $0 OPTIONS

checkBam 0.0.1
author: Colby Chiang (cc2qe@virginia.edu)
description: Check whether a BAM file is truncated. Returns
1 if truncated and 0 if not truncated.

positional arguments:
    bam     BAM file to check

optional arguments:
    -h      Show this message

EOF
}

# Check options passed in.
while getopts "h" OPTION
do
    case $OPTION in
        h)
            usage
            exit 1
            ;;
        ?)
            usage
            exit 1
            ;;
    esac
done

# Check that the arguments are there
if test -z $1
then
    usage
    exit 1
fi

BAMFILE=$1

# get pipe from stdin if -
if [ "$BAMFILE" == "-" ]
then
    BAMFILE="/dev/stdin"
fi

# Do something with the arguments...
OUTPUT=$(( samtools view -H $BAMFILE | head -n 1 1> /dev/null ) 2>&1 )

#echo -e "$OUTPUT"

if [ "$OUTPUT" == "[bam_header_read] EOF marker is absent. The input is probably truncated." ]
then
    echo 1
    exit 1
else
    echo 0
    exit 0
fi

## END SCRIPT
