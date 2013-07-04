#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2013-07-04 15:34 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
gcContentBed.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Generated a BED file of windowed GC content from a fasta file")
    parser.add_argument('-w', '--window', type=int, default=1000, help='window size to calculate GC content over')
    parser.add_argument('-s', '--step', type=int, default=None, help='step size (default: window)')
    parser.add_argument('-b', '--bed', required=False, type=argparse.FileType('w'), default=sys.stdout, help='BED file to output results (default: stdout)')
    parser.add_argument('fasta', nargs='?', type=argparse.FileType('r'), default=None, help='fasta file to read. If \'-\' or absent then defaults to stdin.')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.fasta == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.fasta = sys.stdin

    # if there's no step, make it the window
    if args.step == None:
        args.step = args.window

    # send back the user input
    return args

# primary function
def gcContent(fasta, window, step, bed):
    seqName = None
    seqLength = 0

    prevSeqName = None

    seqcount = 0
    testseq = ''

    for line in fasta:
        line = line.rstrip()

        if len(line) != 0 and line[0] == '>':
            if seqName != None:
                print '%s\t%s' % (seqName, seqLength)
            seqName = line.split()[0].replace('>', '')
            seqLength = 0
        else:
            while seqcount < window:
                getNumber = min(len(line), window-seqcount)
                testseq += line[0:getNumber]
                remainder = line[getNumber:]
                seqcount += getNumber
            print seqcount
            print testseq
            print sum(map(testseq.count, ["G", "C", 'g', 'c']))
            seqcount = 0
            testseq = ''

    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # store into global values
    window = args.window
    step = args.step
    bed = args.bed
    fasta = args.fasta
    
    # call primary function
    gcContent(fasta, window, step, bed)

    # close the input file
    fasta.close()
    bed.close()

# initialize the script
if __name__ == '__main__':
    sys.exit(main())
