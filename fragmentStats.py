#!/usr/bin/env python

import argparse, sys
import pysam
from collections import Counter

# --------------------------------------
# define functions

# get the number of entries in the set
def numEntries(myCounter):
    numEntries = sum(myCounter.values())
    return numEntries

# calculate the mean, given a counter and the
# length of the feature (chromosome or genome)
def mean(myCounter):
    # the number of total entries in the set is the 
    # sum of the occurrences for each value
    numEntries = sum(myCounter.values())
    
    # u holds the mean
    u = float
    mySum = 0.0
    
    for c in myCounter:
        mySum += c * float(myCounter[c])
    u = mySum / numEntries
    return u

def stdev(myCounter):
    # the number of total entries in the set is the 
    # sum of the occurrences for each value
    numEntries = sum(myCounter.values())

    # u holds the mean
    u = mean(myCounter)
    sumVar = 0.0

    # stdev is sqrt(sum((x-u)^2)/#elements)
    for c in myCounter:
        sumVar += myCounter[c] * (c - u)**2
    variance = float(sumVar) / numEntries
    stdev = variance**(0.5)
    return stdev

# calculate quartiles or percentiles from the depth counter.
# percentile(depthCount, 0.5) returns the media
def percentile(myCounter, q):
    # perc is percentile value to return
    perc = float
    #length is the number of bases we're looking at
    numEntries=sum(myCounter.values())
    
    # the number of entries through the set
    # we've gone
    b = 0
    
    # stopping point. Halfway for median, percentile otherwise
    limit = float(numEntries) * q
        
    # a list of the values. Sort it from
    # smallest to largest and start counting to
    # halfway through the # of entries
    valueList = list(myCounter)
    valueList.sort()

    # iterator i
    i = 0
    while b <= limit:
        # if p is exactly halfway, then we must
        # take the mean of the two surrounding.
        if b == limit:
            perc = (valueList[i-1] + valueList[i]) / float(2)
            return perc
    
        myValue = valueList[i]

        # myCounter[myValue] is the number of entries
        # that have myValue
        b += myCounter[myValue]
        i += 1
    perc = valueList[i-1]
    return perc

# median is simply 50th percentile
def median(myCounter):
    return percentile(myCounter, 0.5)

# Driver function
def fragmentStats(bamfile):
    if bamfile.name == '<stdin>':
        myBam = pysam.Samfile('-', 'rb')
    else:
        myBam = pysam.Samfile(bamfile.name, 'rb')

    # counters to hold the depth counts over the chromosome
    # and the whole genome
    mappedFrags = Counter()

    for samRec in myBam:

#        if samRec.is_proper_pair and samRec.is_read1:
        if samRec.is_proper_pair and not samRec.is_reverse:
            mappedFrags[abs(samRec.tlen)] += 1

    print 'median ', median(mappedFrags)
    print 'mean ', mean(mappedFrags)
    print 'stdev ', stdev(mappedFrags)
    print 'numRecs ', numEntries(mappedFrags)
    

    return

# --------------------------------------
# argument parsing

def main():
    parser = argparse.ArgumentParser(description="Calculated fragment length statistics from a coverage SAM file.")

    parser.add_argument('-i', '--input', required=True, type=argparse.FileType('r'), help='SAM file (\'-\' for stdin)')
    
    # parse the arguments
    args = parser.parse_args()

    # store into global values
    bamfile = args.input
    
    # call the driver function
    fragmentStats(bamfile)

    # close the file
    bamfile.close()

# initialized the script
if __name__ == '__main__':
    sys.exit(main())
