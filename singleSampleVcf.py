#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2013-11-22 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
singleSampleVcf.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Basic python script template")
    parser.add_argument('vcf', nargs='?', type=argparse.FileType('r'), default=None, help='vcf file to read. If \'-\' or absent then defaults to stdin.')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.vcf == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.vcf = sys.stdin

    # send back the user input
    return args

# primary function
def singleSampleVcf(vcf):
    for line in vcf:
        if line[0:2] == '##':
            print line.rstrip()
            continue
        elif line[0:6] == '#CHROM':
            v = line.rstrip().split('\t')
            print '\t'.join(v[:10])
            continue

        v = line.rstrip().split('\t')
        (chrom, pos, varid, ref, alt, qual, filt, info, format, sampleInfo) = v
        alt_split = alt.split(',')
        if len(alt_split) == 1:
            print line.rstrip()
        else:
            sampleInfo_split = sampleInfo.split(':')
            gt = sampleInfo_split[0]
            if gt == '.':
                continue
            gt_split = gt.split('/')
            gt_split = map(int, gt_split)

            # print "OLD", chrom, pos, ref, alt_split, gt_split

            new_alt_split = list()
            new_gt_split = list()
            for g in gt_split:
                alt_allele = alt_split[g-1]
                if g != 0 and alt_allele not in new_alt_split:
                    new_alt_split.append(alt_allele)
                if g != 0:
                    new_gt = new_alt_split.index(alt_allele) + 1
                    new_gt_split.append(new_gt)
                elif g == 0:
                    new_gt_split.append(0)

            new_alt = ','.join(new_alt_split)
            new_gt = '/'.join(map(str, new_gt_split))
            new_sampleInfo = new_gt + ':' + ':'.join(sampleInfo_split[1:])
            print '\t'.join(map(str, (chrom, pos, varid, ref, new_alt, qual, filt, info, format, new_sampleInfo)))

            
    
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # store into global values
    vcf = args.vcf
    
    # call primary function
    singleSampleVcf(vcf)

    # close the input file
    vcf.close()

# initialize the script
if __name__ == '__main__':
    sys.exit(main())
