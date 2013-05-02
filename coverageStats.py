#!/usr/bin/env python

import argparse, sys
from collections import Counter

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
    print '\t'.join(('chr', 'length', 'basesCovered', 'portionCovered', 'meanDepth', 'stdevDepth', 'medianDepth', 'q1', 'q3'))

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

                basesCov = basesCovered(chromDepth)
                fracCov = float(basesCov) / chromLength
                print '\t'.join(map(str,(prevChrom, chromLength, basesCov, '%.3g' % fracCov, '%.3g' % mean(chromDepth, chromLength), '%.3g' % stdev(chromDepth, chromLength), int(median(chromDepth)), int(percentile(chromDepth, 0.25)), int(percentile(chromDepth, 0.75)) )))
                
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
    basesCov = basesCovered(chromDepth)
    fracCov = float(basesCov) / chromLength
    print '\t'.join(map(str,(prevChrom, chromLength, basesCov, '%.3g' % fracCov, '%.3g' % mean(chromDepth, chromLength), '%.3g' % stdev(chromDepth, chromLength), int(median(chromDepth)), int(percentile(chromDepth, 0.25)), int(percentile(chromDepth, 0.75)) )))

    # print the whole genome coverage statistics    
    genomeLength += end
    basesCov = basesCovered(genomeDepth)
    fracCov = float(basesCov) / genomeLength
    print '\t'.join(map(str,('total', genomeLength, basesCov, '%.3g' % fracCov, '%.3g' % mean(genomeDepth, genomeLength), '%.3g' % stdev(genomeDepth, genomeLength), int(median(genomeDepth)), int(percentile(genomeDepth, 0.25)), int(percentile(genomeDepth, 0.75)) )))
    return

# --------------------------------------
# argument parsing

def main():
    parser = argparse.ArgumentParser(description="Calculated read depth statistics from a coverage BED file.")

    parser.add_argument('-i', '--input', required=True, type=argparse.FileType('r'), help='Tab delimited BED file that reports regions of zero coverage (\'-\' for stdin)')
    
    # parse the arguments
    args = parser.parse_args()

    # store into global values
    bedfile = args.input
    
    coverageStats(bedfile)
    bedfile.close()

# initialized the script
if __name__ == '__main__':
    sys.exit(main())
