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
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), default=sys.stdin, help='Fasta file to read (default: stdin)')
    
    # parse the arguments
    args = parser.parse_args()

    # store into global values
    fasta = args.input

    # call the driver function
    getSeqLengths(fasta)

    # close the file
    fasta.close()

# initialize the script
if __name__ == '__main__':
    sys.exit(main())
