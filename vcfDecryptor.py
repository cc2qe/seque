#!/usr/bin/env python

import argparse, sys, re
from argparse import RawTextHelpFormatter

__author__ = "Author (email@site.com)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2013-05-09 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, add_help=False, description="\
pythonTemplate.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Basic python script template")
    #    parser.add_argument('-a', '--argA', metavar='argA', type=str, required=True, help='description of argument')
    # parser.add_argument('-b', '--argB', metavar='argB', required=False, help='description of argument B')
    parser.add_argument('-h', '--header', required=False, action='store_true', help='print header')
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
def myFunction(vcf, header):
    ann_set = set()
    header_done = False
    for l in vcf:
        line = l.rstrip()
        if line[0] != '#':
            if not header_done:
                headOut = ['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter']
                for ann in ann_set:
                    headOut.append(ann)
                print '\t'.join(headOut)
                header_done = True
            (chrom, pos, varId, ref, alt, qual, filt, info_string, format) = line.split('\t')[:9]
            pos = int(pos)

            outList = [chrom, pos, varId, ref, alt, qual, filt]

            info_list = info_string.split(';')
            info_dict = dict()
            for i in info_list:
                if len(i.split('=')) > 1:
                    (ann, value) = i.split('=')
                    info_dict[ann] = value
            for ann in ann_set:
                if ann in info_dict:
                    outList.append(info_dict[ann])
                else:
                    outList.append('na')

            print '\t'.join(map(str, outList))


        else:
            if line[:6]=='##INFO':
                info_field = re.sub('^##INFO=<ID=', '', line).split(',')[0]
                ann_set.add(info_field)
                
            if header:
                print line
    
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # store into global values
    vcf = args.vcf
    header = args.header
    
    # call primary function
    myFunction(vcf, header)

    # close the input file
    vcf.close()

# initialize the script
if __name__ == '__main__':
    sys.exit(main())
