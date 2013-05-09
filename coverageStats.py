#!/usr/bin/env python

import argparse, sys
from collections import Counter
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2013-05-09 14:31 $"


# --------------------------------------
# define functions

# calculate the mean, given a counter and the
# length of the feature (chromosome or genome)
def mean(depthCount, length):
    # u holds the mean
    u = float
    sumCov = 0.0
    
    for d in depthCount:
        sumCov += d * float(depthCount[d])
    u = sumCov / length
    return u

def stdev(depthCount, length):
    # u holds the mean
    u = mean(depthCount, length)
    sumVar = 0.0

    # stdev is sqrt(sum((x-u)^2)/#elements)
    for d in depthCount:
        sumVar += (d - u)**2
    variance = float(sumVar) / length
    stdev = variance**(0.5)
    return stdev

# calculate quartiles or percentiles from the depth counter.
# percentile(depthCount, 0.5) returns the media
def percentile(depthCount, q):
    # perc is percentile value to return
    perc = float
    #length is the number of bases we're looking at
    length=sum(depthCount.values())
    
    # the number of bases through the region
    # we've gone
    b = 0
    
    # stopping point. Halfway for median
    limit = float(length) * q
        
    # a list of the depth numbers. Sort it from
    # smallest to largest and start counting to
    # halfway through the # of bases
    depthList = list(depthCount)
    depthList.sort()

    # iterator i
    i = 0
    while b <= limit:
        # if p is exactly halfway, then we must
        # take the mean of the two surrounding.
        if b == limit:
            perc = (depthList[i-1] + depthList[i]) / float(2)
            return perc
    
        # myDepth is the coverage depth
        myDepth = depthList[i]

        # depthCount[myDepth] is the number of bases
        # that have depth of myDepth
        b += depthCount[myDepth]
        i += 1
    perc = depthList[i-1]
    return perc

# median is simply 50th percentile
def median(depthCount):
    return percentile(depthCount, 0.5)

# number of bases where coverage is greater than zero
def basesCovered(depthCount):
    length = sum(depthCount.values())
    bCov = length - depthCount[0]
    return bCov

# Driver function
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
    print '\t'.join(('chr', 'length', 'basesCovered', 'portionCovered', 'mean', 'stdev', 'median', 'q1', 'q3', '2.5%', '97.5%', 'min', 'max', 'range'))

    for line in bedfile:
        v = line.rstrip().split('\t')
        chrom = v[0]

        # if you run into a new chrom, print the
        # stats on the previous chrom
        if chrom != prevChrom:
            # first make sure you're not looking at the
            # first chrom
            if prevChrom != None:
                # calc bases and fraction of chrom covered.
                basesCov = basesCovered(chromDepth)
                fracCov = float(basesCov) / chromLength

                # print tabulated data for previous chrom
                print '\t'.join(map(str,(prevChrom, chromLength, basesCov, fracCov, mean(chromDepth, chromLength), stdev(chromDepth, chromLength), int(median(chromDepth)), int(percentile(chromDepth, 0.25)), int(percentile(chromDepth, 0.75)), int(percentile(chromDepth, 0.025)), int(percentile(chromDepth, 0.975)), min(list(chromDepth)), max(list(chromDepth)), max(list(chromDepth)) - min(list(chromDepth)) )))
                # start a new counter
                chromDepth=Counter()

                # set new chromLength to zero
                chromLength = 0
                
        # cast the start end, and depth
        (start, end) = map(int, v[1:3])
        depth = float(v[3])
        span = (end-start)

        # add span to depth counts
        chromDepth[depth] += span
        genomeDepth[depth] += span

        # add span to the chrom and genome sizes
        chromLength += span
        genomeLength += span

        # store the previous line's data
        prevChrom = chrom

    # calc bases and fraction of chrom covered
    basesCov = basesCovered(chromDepth)
    fracCov = float(basesCov) / chromLength

    # print tabulated data for last chrom
    print '\t'.join(map(str,(prevChrom, chromLength, basesCov, fracCov, mean(chromDepth, chromLength), stdev(chromDepth, chromLength), int(median(chromDepth)), int(percentile(chromDepth, 0.25)), int(percentile(chromDepth, 0.75)), int(percentile(chromDepth, 0.025)), int(percentile(chromDepth, 0.975)), min(list(chromDepth)), max(list(chromDepth)), max(list(chromDepth)) - min(list(chromDepth)) )))
                
    # print the whole genome coverage statistics    
    basesCov = basesCovered(genomeDepth)
    fracCov = float(basesCov) / genomeLength

    print '\t'.join(map(str,('total', genomeLength, basesCov, fracCov, mean(genomeDepth, genomeLength), stdev(genomeDepth, genomeLength), int(median(genomeDepth)), int(percentile(genomeDepth, 0.25)), int(percentile(genomeDepth, 0.75)), int(percentile(genomeDepth, 0.025)), int(percentile(genomeDepth, 0.975)), min(list(genomeDepth)), max(list(genomeDepth)), max(list(genomeDepth)) - min(list(genomeDepth)) )))
    return

# --------------------------------------
# argument parsing

def main():
    parser = argparse.ArgumentParser(description="coverageStats.py\n\
author: Colby Chiang (cc2qe@virginia.edu)\n\
version: 0.0.1\n\
description: Calculate read depth statistics from a coverage BED file.", formatter_class=RawTextHelpFormatter)

    parser.add_argument('bedfile', nargs='?', type=argparse.FileType('r'), default=None, help='Tab delimited BED file that reports regions of zero coverage (default: stdin)')
    
    # parse the arguments
    args = parser.parse_args()

    # store into global values
    bedfile = args.bedfile

    # if no bedfile, check if part of pipe and if so, read stdin.
    if bedfile == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            bedfile = sys.stdin
    
    coverageStats(bedfile)
    bedfile.close()

# initialized the script
if __name__ == '__main__':
    sys.exit(main())
