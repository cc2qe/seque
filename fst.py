#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter
import vcf as pyvcf

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2013-08-01 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
fst.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Calculate Fst")
    parser.add_argument('-s', '--samples', type=argparse.FileType('r'), required=True, help='sample file')
    parser.add_argument('vcfFile', metavar='vcf', nargs='?', type=argparse.FileType('r'), default=None, help='tab-delimited variant file to read. If \'-\' or absent then defaults to stdin.')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.vcfFile == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.vcfFile = sys.stdin

    # send back the user input
    return args

# primary function
def getFst(sampFile, vcfFile):
    
    # for line in varFile:
    #    print line.rstrip()
    vcf_reader = pyvcf.Reader(open(vcfFile.name), 'r', compressed=True)
    #vcf_reader = pyvcf.Reader('-', 'r')

    for record in vcf_reader:
        print record.ID
        for sample in record.samples:
            print sample['GT']
    
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # store into global values
    sampFile = args.samples
    vcfFile = args.vcfFile
    
    # call primary function
    getFst(sampFile, vcfFile)

    # close the input file
    vcfFile.close()
    sampFile.close()

# initialize the script
if __name__ == '__main__':
    sys.exit(main())
