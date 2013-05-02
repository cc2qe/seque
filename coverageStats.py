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
        sumCov += d * float(depthCount[d])
    mean = sumCov / length
    return mean

def median(depthCount):
    median = float
    length=sum(depthCount.values())
    
    # the number of bases through the region
    # we've gone with median
    p = 0
    
    #
    halfway = float(length) / 2
        
    # a list of the depth numbers. Sort it from
    # smallest to largest and start counting to
    # halfway through the # of bases
    depthList = list(depthCount)
    depthList.sort()

    # iterator i
    i = 0
    while p <= halfway:
        # if p is exactly halfway, then we must
        # take the mean of the two surrounding.
        if p == halfway:
            median = (depthList[i-1] + depthList[i]) / float(2)
            return median
    
        # myDepth is the coverage depth
        myDepth = depthList[i]

        # depthCount[myDepth] is the number of bases
        # that have depth of myDepth
        p += depthCount[myDepth]
        i += 1
    median = depthList[i-1]
    return median

def coverageStats(bedfile):
    
    # vars for storing the previous line's info
    prevChrom = None
    prevEnd = 0

    # counters to hold the depth counts over the chromosome
    # and the whole genome
    chromDepth = Counter()
    genomeDepth = Counter()

    # initialize the chrom and genome lengths to zero
    chromLength = 0
    genomeLength = 0

    # declare global vars for each bed line
    chrom = str
    start = int
    end = int
    depth = float
    
    # print header
    print '\t'.join(('chr', 'length', 'meanDepth', 'medianDepth'))

    for line in bedfile:
        v = line.rstrip().split('\t')
        chrom = v[0]

        # if you run into a new chrom, print the
        # stats on the previous chrom
        if chrom != prevChrom:
            # first make sure you're not looking at the
            # first chrom
            if prevChrom != None:
                # set the chrom length and iterate the genome length
                # (assumes BED file has all positions in each chrom
                # even if it's zero coverage
                chromLength = prevEnd
                genomeLength += prevEnd
                print '\t'.join(map(str,(prevChrom, chromLength, mean(chromDepth, chromLength), median(chromDepth))))

                # start a new counter
                chromDepth=Counter()

        # cast the start end, and depth
        (start, end) = map(int, v[1:3])
        depth = float(v[3])
        span = (end-start)

        # add span to depth counts
        chromDepth[depth] += span
        genomeDepth[depth] += span
        
        # store the previous line's data
        prevChrom = chrom
        prevEnd = end


    # print the last chromosome
    chromLength = end
    print '\t'.join(map(str,(chrom, chromLength, mean(chromDepth, chromLength), median(chromDepth))))



    # print the whole genome coverage statistics    
    genomeLength += end
    print '\t'.join(map(str,('total', genomeLength, mean(genomeDepth, genomeLength), median(genomeDepth))))
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
