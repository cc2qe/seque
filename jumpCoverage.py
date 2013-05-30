#!/usr/bin/env python


import sys
import argparse
import pdb
import pysam


def find_insert(pos, regions):
    for region in regions:
        if pos >= region.start and pos <= region.end:
            return region


class Region():
    def __init__(self, chr, start, end, coverage):
        self.chr = chr
        self.start = start
        self.end = end
        self.cov = coverage

    def split(self, region):
        left = Region(self.chr, self.start, region.start, self.cov)
        mid = Region(self.chr, region.start, region.end, self.cov + region.cov)
        right = Region(self.chr, region.end, self.end, self.cov)

        # Don't need to check for start == end because that's just 1 bp region
        #return [left, mid, right]
        splits = [left, mid, right]
        return [r for r in splits if (r.end - r.start) > 0]

    def split_left(self, region):
        left = Region(self.chr, self.start, region.start, self.cov)
        right = Region(self.chr, region.start, self.end, self.cov + region.cov)
        #return [left, right]
        return [r for r in [left, right] if (r.end - r.start) > 0]

    def split_right(self, region):
        left = Region(self.chr, self.start, region.end, self.cov + region.cov)
        right = Region(self.chr, region.end, self.end, self.cov)
        #return [left, right]
        return [r for r in [left, right] if (r.end - r.start) > 0]

    def __repr__(self):
        return "%s\t%s\t%s\t%s" % (self.chr, self.start, self.end, self.cov)


class BedChrom():
    def __init__(self, region):
        self.regions = [region]

    def insert_region(self, new_region):
        left_insert = find_insert(new_region.start, self.regions)
        right_insert = find_insert(new_region.end, self.regions)

        # Sometimes in mitochondrial DNA and the pseudo chromosomes
        # (chrUn, chrN_random) reads will map closer to the end of the
        # pseudo chromosome than the length of read. In these cases,
        # find_insert will fail, since the end position is off the end of the
        # chromosome, so we use the end region.
        if not right_insert:
            right_insert = self.regions[-1]

        if left_insert == right_insert:
            i = self.regions.index(left_insert)
            self.regions = self.regions[:i] + left_insert.split(new_region) + self.regions[i + 1:]
        else:
            j = self.regions.index(left_insert)
            k = self.regions.index(right_insert)

            for x in xrange(j + 1, k):
                self.regions[x].cov += new_region.cov

            self.regions = (self.regions[:j] + left_insert.split_left(new_region) +
                self.regions[j + 1:k] + right_insert.split_right(new_region) +
                self.regions[k + 1:])

    def write_out(self, outfile, pos):
        for region in self.regions:
            if region.start < pos:
                outfile.write("%s\n" % str(region))
            else:
                i = self.regions.index(region)
                self.regions = self.regions[i:]
                break


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate a BED graph of "
        "jumping library coverage.")
    parser.add_argument('samfile', nargs='?', type=pysam.Samfile,
        default=sys.stdin, help="Input samfile. Defaults to stdin.")
    parser.add_argument('bedgraph', nargs='?', type=argparse.FileType('w'),
        default=sys.stdout, help="Output bedgraph file. Defaults to stdout.")
    #parser.add_argument('genome', type=argparse.FileType('r'),
    #    help="Tab delimited file of chromosome sizes in genome\n"
    #    "<chromName><TAB><chromSize>")

    args = parser.parse_args()

    samfile = args.samfile
    bedgraph = args.bedgraph

    genome = {}
    for contig in samfile.header['SQ']:
        genome[contig['SN']] = contig['LN']

    #genome = []
    #for line in args.genome:
    #    genome.append(line.strip().split('\t'))
    #genome = dict(genome)

    curr_chr = ''
    bed_chrom = None

#    pdb.set_trace()
    for read in samfile:
        # Skip unmapped reads
        if 0 <= read.tid < samfile.nreferences:
            chrom_name = samfile.getrname(read.tid)
        else:
            continue

        # If new chromosome
        if chrom_name != curr_chr:
            # Print everything remaining in previous chromosome
            if bed_chrom:
                bed_chrom.write_out(bedgraph, sys.maxint)
            # Begin new chromosome
            # print "Beginning coverage calculation over %s" % chrom_name
            curr_chr = chrom_name
            chr_size = int(genome[curr_chr])
            bed_chrom = BedChrom(Region(curr_chr, 0, chr_size, 0))

        # Add it
        if read.is_proper_pair and read.tlen > 0:
            # Insert region into current list
            # Note: pysam and bed are both zero based, samtools is one based
            # Region is half-open interval
            bed_chrom.insert_region(Region(curr_chr, read.pos, read.mpos + read.qlen, 1))

            # Print regions up to start of this read
            bed_chrom.write_out(bedgraph, read.pos)
