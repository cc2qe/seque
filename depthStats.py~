#!/usr/bin/env python

import argparse, sys

# --------------------------------------
# define functions

def myFunction(bedfile):
    
    genomeSum = 0.0
    genomeNumBases = 0

    currentChrom = None
    chromSum = 0.0
    chromNumBases = 0
    
    # print header
    print '\t'.join(('chr', 'length', 'mean depth'))

    for line in bedfile:
        v = line.rstrip().split('\t')

        chrom = v[0]
        (start, end) = map(int, v[1:3])
        depth = float(v[3])
        
        span = (end-start)
        genomeSum += depth * span
        chromSum += depth * span
        genomeNumBases += span
        chromNumBases += span

        # print chrom, depth
        if chrom != currentChrom:
            if currentChrom != None:
                print '\t'.join(map(str,(currentChrom, chromNumBases, (chromSum/chromNumBases))))
            currentChrom = chrom
            chromSum = 0.0
            chromNumBases = 0

    # print the last chromosome
    print '\t'.join(map(str,(currentChrom, chromNumBases, (chromSum/chromNumBases))))

    # print the overall stats
    print '\t'.join(map(str,('Total', genomeNumBases, (genomeSum/genomeNumBases))))
    return

# --------------------------------------
# argument parsing

def main():
    parser = argparse.ArgumentParser(description="Calculated read depth statistics from a coverage BED file.")

    parser.add_argument('bedfile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='Tab delimited BED file that reports regions of zero coverage (default: stdin)')
    
    # parse the arguments
    args = parser.parse_args()

    # store into global values
    bedfile = args.bedfile
    
    myFunction(bedfile)
    bedfile.close()

# initialized the script
if __name__ == '__main__':
    sys.exit(main())
