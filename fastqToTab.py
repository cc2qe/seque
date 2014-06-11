#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cc2qe@virginia.edu"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2013-11-14 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
fastqToTab.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: convert a fastq file into a 4 column tab delimited file")
    parser.add_argument('-i', '--interleaved', required=False, action='store_true', help='fastq file is interleaved, with pairs consecutive (outputs 8 line tab delimited file')
    parser.add_argument('fastq', nargs='?', type=argparse.FileType('r'), default=None, help='fastq file to read. If \'-\' or absent then defaults to stdin.')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.fastq == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.fastq = sys.stdin

    # send back the user input
    return args

# primary function
def fastqToTab(fastq, interleaved):
    counter = 1

    # if interleaved then output batches of 8 columns with both paired-end reads
    if interleaved:
        interval = 8
    else:
        interval = 4

    for line in fastq:
        if counter % interval == 0:
            sys.stdout.write(line)
        else:
            sys.stdout.write(line.rstrip() + "\t")

        counter = counter + 1
        
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # store into global values
    fastq = args.fastq
    interleaved = args.interleaved
    
    # call primary function
    fastqToTab(fastq, interleaved)

    # close the input file
    fastq.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
