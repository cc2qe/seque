#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter

__author__ = "Author (email@site.com)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2013-05-09 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
covDepth.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Calculate the coverage depth over the GRCh37 human genome")
    parser.add_argument('-n', '--numreads', type=int, required=True, help='number of reads')
    parser.add_argument('-l', '--readlength', type=int, required=True, help='length of each read')

    # parse the arguments
    args = parser.parse_args()

    # send back the user input
    return args

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    genome_size = 3101804739
    
    cov = args.numreads * args.readlength / float(genome_size)

    print cov
    


# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
