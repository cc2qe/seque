#!/usr/bin/env python

import argparse, sys

# --------------------------------------
# define functions

def getSeqLengths(fasta):
    seqName = None
    seqLength = 0

    prevSeqName = None
    
    for line in fasta:
        line = line.rstrip()
        if len(line) != 0 and line[0] == '>':
            if seqName != None:
                print '%s\t%s' % (seqName, seqLength)
            seqName = line.split()[0].replace('>', '')
            seqLength = 0
        else:
            seqLength += len(line)
    print '%s\t%s' % (seqName, seqLength)
    return

# --------------------------------------
# argument parsing

def main():
    parser = argparse.ArgumentParser(description="Read a fasta file and output a tab delimited file of sequence names and lengths.")
    parser.add_argument('fasta', nargs='?', type=argparse.FileType('r'), default=None, help='Fasta file to read (default: stdin)')
    
    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.fasta == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.fasta = sys.stdin

    # store into global values
    fasta = args.fasta

    # call the driver function
    getSeqLengths(fasta)

    # close the file
    fasta.close()

# initialize the script
if __name__ == '__main__':
    sys.exit(main())
