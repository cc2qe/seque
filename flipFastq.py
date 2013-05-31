#!/usr/bin/env python

import sys, argparse
from string import *

# gets the reverse complement of a DNA strand
def revcomp(dna):
    """ reverse complement of a DNA sequence """
    
    comp = dna.translate(maketrans("AGCTagct", "TCGAtcga"))
    lcomp = list(comp)
    lcomp.reverse()

    return ''.join(lcomp)


# argument parsing
# -----------------------

def main():
    parser = argparse.ArgumentParser(description='Flip fastq reads and their base qualities to output the reverse complement reads. Outputs the flipped fastq file to stdout. Useful for changing outward facing reads to inward facing or vice versa')
    parser.add_argument('input', nargs='?', type=argparse.FileType('r'), default=None, help='Input fastq file. (default: stdin)')

    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.input == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.input = sys.stdin

    f = args.input

    i = 1
    for l in f:
	l = l.rstrip()

        # revcomp the sequence
	if i % 4 == 2:
            print revcomp(l)
                
        # reverse the base quality string
	elif i % 4 == 0:
            print l[::-1]
        
        # print the read names lines as is
	else:
            print l
	i += 1

# initialize the script
if __name__ == '__main__':
    sys.exit(main())
