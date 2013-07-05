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
    parser.add_argument('-w', '--window', type=int, required=True, help='window size to calculate GC content over')
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
    # fasta should start with a >chrom line, so we'll grab the first one
    line = fasta.readline().rstrip()
    if line[0] == '>':
        chrom = line[1:]
    else:
        sys.stderr.write("Error: fasta does not begin with a >chrom line\n")
    nextChrom = chrom

    seqLen = 0
    seq = ''
    remainder = ''
    startPos = 0

    while 1:
        while seqLen < window and len(remainder) != 0:
            # before reading the next line, burn through the remainder
            getNumber = min(len(remainder), window-seqLen)
            seq += remainder[0:getNumber]
            remainder = remainder[getNumber:]
            seqLen += getNumber

        while seqLen < window:
            line = fasta.readline().rstrip()
            if line == '':
                break
            
            if line[0] == '>':
                nextChrom = line[1:]
                break

            line = remainder + line
            #print 'line with rem', line

            getNumber = min(len(line), window-seqLen)
            seq += line[0:getNumber]
            remainder = line[getNumber:]
            seqLen += getNumber

        #print 'seq', seq
        #print 'seqlen', seqLen
        #print 'rem', remainder
        gcCount = sum(map(seq.count, ['G', 'C', 'g', 'c']))
        #print gcCount

        gcFrac = gcCount/float(seqLen)
        # maybe this should be min startPos+step or end of chrom            
        bed.write('\t'.join(map(str, [chrom, startPos, startPos+step, gcFrac])) + '\n')
    
        startPos += step
        # careful of this at the end of a chromosome
        if chrom == nextChrom:
            remainder = seq[step:] + remainder # i think it's still possible to have a remainder at this point
        chrom = nextChrom
        #print 'rem at bottom', remainder
        seqLen = 0
        seq = ''

        if line == '':
            break
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
