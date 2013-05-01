#!/usr/bin/env python

import argparse, sys
from collections import Counter

# --------------------------------------
# define functions

def myFunction(bedfile):
    
    currentChrom = None

    chromDepth = Counter()
    genomeDepth = Counter()
    
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
                chromSize = 0
                for d in chromDepth:
                    sumCov += d * chromDepth[d]
                    chromSize += chromDepth[d]
                print '\t'.join(map(str,(currentChrom, chromSize, (sumCov/chromSize))))
            currentChrom = chrom
            chromDepth=Counter()

    # print the last chromosome
    sumCov = 0
    chromSize = 0
    for d in chromDepth:
        sumCov += d * chromDepth[d]
        chromSize += chromDepth[d]
    print '\t'.join(map(str,(currentChrom, chromSize, (sumCov/chromSize))))

    # print the overall stats
    # print '\t'.join(map(str,('Total', genomeNumBases, (genomeSum/genomeNumBases))))

    # print genomeDepth
    
    totalBases = 0
    spanSum = 0
    for d in genomeDepth:
        totalBases += (d * genomeDepth[d])
        spanSum += genomeDepth[d]

    print totalBases, spanSum
    print 'by counts, mean is %s' % (totalBases/spanSum)
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
    
    myFunction(bedfile)
    bedfile.close()

# initialized the script
if __name__ == '__main__':
    sys.exit(main())
