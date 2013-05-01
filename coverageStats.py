#!/usr/bin/env python

import argparse, sys
from collections import Counter

# --------------------------------------
# define functions

# calculate the mean, given a counter and the
# length of the feature (chromosome or genome)
def mean(depthCount, length):
    mean = float
    sumCov = 0.0
    
    for d in depthCount:
        sumCov += d * float(chromDepth[d])
    mean = sumCov / length
    
    return mean

def coverageStats(bedfile):
    
    currentChrom = None

    chromDepth = Counter()
    genomeDepth = Counter()

    chromSize = 0
    genomeSize = 0

    # declare global vars for each bed line
    chrom = str
    start = int
    end = int
    depth = float
    
    # print header
    print '\t'.join(('chr', 'length', 'mean depth'))

    for line in bedfile:
        v = line.rstrip().split('\t')

        chrom = v[0]
        (start, end) = map(int, v[1:3])
        depth = float(v[3])
        
        span = (end-start)

        # add span to depth counts
        chromDepth[depth] += span
        genomeDepth[depth] += span

        # print chrom, depth
        if chrom != currentChrom:
            if currentChrom != None:
                
                
                sumCov = 0
                chromSize = end
                genomeSize += end
                for d in chromDepth:
                    sumCov += d * chromDepth[d]
                print '\t'.join(map(str,(currentChrom, chromSize, (sumCov/chromSize))))
            currentChrom = chrom
            chromDepth=Counter()

    # print the last chromosome
    sumCov = 0
    chromSize = 0
    for d in chromDepth:
        sumCov += d * chromDepth[d]
        chromSize = end
    print '\t'.join(map(str,(currentChrom, chromSize, (sumCov/chromSize))))

    # print the whole genome coverage statistics    
    sumCov = 0
    for d in genomeDepth:
        sumCov += (d * genomeDepth[d])
        genomeSize += genomeDepth[d]
    print '\t'.join(map(str,('total', genomeSize, (sumCov/genomeSize))))
    return

# --------------------------------------
# argument parsing

def main():
    parser = argparse.ArgumentParser(description="Calculated read depth statistics from a coverage BED file.")

    parser.add_argument('-i', '--input', type=argparse.FileType('r'), help='Tab delimited BED file that reports regions of zero coverage (\'-\' for stdin)')
    
    # parse the arguments
    args = parser.parse_args()

    # store into global values
    bedfile = args.input
    
    coverageStats(bedfile)
    bedfile.close()

# initialized the script
if __name__ == '__main__':
    sys.exit(main())
