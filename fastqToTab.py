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
def fastqToTab(fastq):
    counter = 1
    for line in fastq:

        if counter % 4 == 0:
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
    
    # call primary function
    fastqToTab(fastq)

    # close the input file
    fastq.close()

# initialize the script
if __name__ == '__main__':
    sys.exit(main())
