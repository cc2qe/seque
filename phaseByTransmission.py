#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2014-01-29 21:35 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
phaseByTransmission.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Phase a VCF file based on a pedigree")
    parser.add_argument('-r', '--retain', action='store_true', help='retain original phasing (default: overwrite)')
    parser.add_argument('-f', '--fatherFirst', action='store_true', help='outputs phasing as father|mother (default is mother|father)')
    parser.add_argument('-p', '--pedigree', type=argparse.FileType('r'), required=True, help='PED format file (http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped)') 
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

# vcf class
class Vcf(object):
    def __init__(self):
        self.header = []
        self.col_headers = []
        self.sample_list = []
        return
    def add_header_info(self, header_line):
        self.header.append(header_line)
        if header_line[0:2] != '##' and header_line[0] == '#':
            v = header_line.split('\t')
            self.col_headers = v
            self.sample_list = self.col_headers[9:]
        return
    def print_header(self):
        for l in self.header:
            print l
        return

# class of each variant within the VCF file
class Variant(object):
    def __init__(self, myVcf, vcf_line):
        # set the Vcf
        self.vcf = myVcf
        self.sample_list = myVcf.sample_list
        self.num_samples = len(self.sample_list)
        
        # parse the vcf line
        self.vcf_line = vcf_line
        v = vcf_line.split('\t')
        self.chrom = v[0]
        self.pos = int(v[1])
        self.id = v[2]
        self.ref = v[3]
        self.alt = v[4].split(',')
        self.qual = float(v[5])
        self.filter = v[6]
        self.info = v[7]
        self.format = v[8].split(':')
        
        # store the gt field for each sample
        #self.sample_info = [None] * self.num_samples
        self.sample_dict = dict()
        for i in xrange(self.num_samples):
            self.sample_dict[self.sample_list[i]] = Genotype(self.sample_list[i], v[9+i], self)

    def sample_info(self, samp):
        return self.sample_dict[samp]

    # print the variant line
    def get_var_string(self):
        var_fields = [self.chrom, self.pos, self.id, self.ref, ','.join(self.alt), self.qual, self.filter, self.info, ':'.join(self.format)]
        for s in self.sample_list:
            var_fields.append(self.sample_info(s).get_field_string())
        return '\t'.join(map(str, var_fields))

# each sample has a genotype at each variant, which includes
# its genotype as well as allele depth stats, et al.
class Genotype(object):
    def __init__(self, sample_id, gt_field, myVar):
        self.sample_id = sample_id
        self.myVar = myVar

        self.format = myVar.format
        self.gt_list = gt_field.split(':')
        
        # mark as phased or unphased and store alleles
        if '|' in self.gt_list[0]:
            self.phased = True
            self.alleles = self.gt_list[0].split('|')
        else:
            self.phased = False
            self.alleles = self.gt_list[0].split('/')

    def is_phased(self):
        return self.phased

    def get_alleles(self):
        return self.alleles

    def get_gt_string(self):
        if self.is_phased():
            delim = '|'
        else:
            delim = '/'
        return delim.join(self.alleles)

    def get_field_string(self):
        return ':'.join([self.get_gt_string()] + self.gt_list[1:])

# a Pedigree is a dictionary of ped entries where the key is the sample
class Ped_entry(object):
    def __init__(self, v):
        self.fam = v[0]
        self.sample = v[1]
        self.pat = v[2]
        self.mat = v[3]
        self.sex = v[4]
        self.phen = v[5]
        if len(v) > 6:
            self.other = v[7:]

class Pedigree(object):
    # init the Pedigree
    def __init__(self, pedFile):
        # make a list of pedigree objects for each of the samples
        self.ped = dict()
        self.sample_list = []
        for line in pedFile:
            p = Ped_entry(line.rstrip().split('\t'))
            self.ped[p.sample] = p
            self.sample_list.append(p.sample)

    def get_sample(self, sample):
        return self.ped[sample]
        
# primary function
def phase_by_transmission(myVar, pedigree, fatherFirst, retain):
    for samp_id in pedigree.sample_list:

        samp_ped = pedigree.get_sample(samp_id)
        mom_id = samp_ped.mat
        dad_id = samp_ped.pat

        # set samp
        if samp_id in myVar.sample_list:
            samp = myVar.sample_info(samp_id)

        # skip if sample is already phased and retain flag is set
        if retain and samp.is_phased():
            continue

        # store a and b as the alleles in the sample
        if len(samp.alleles) == 2 and '.' not in samp.alleles:
            a = samp.alleles[0]
            b = samp.alleles[1]
        else:
            continue

        # homozygous vars are always in phase
        if a == b:
            samp.phased = True
            continue

        # check that the mom and dad are in the VCF
        if mom_id in myVar.sample_list and dad_id in myVar.sample_list:
            mom = myVar.sample_info(mom_id)
            dad = myVar.sample_info(dad_id)
        else:
            continue

        # skip if parent has nocall
        if '.' in mom.alleles + dad.alleles:
            continue

        # check for mendelian transmission
        #if (a not in (mom.alleles + dad.alleles)) or (b not in (mom.alleles + dad.alleles)):
            #print "mendelian error: %s: %s, %s: %s, %s: %s" % (samp_id, samp.get_gt_string(), mom_id, mom.get_gt_string(), dad_id, dad.get_gt_string())


        if a in mom.alleles and a not in dad.alleles and b in dad.alleles:
            if fatherFirst:
                samp.alleles = [b, a]
            else:
                samp.alleles = [a, b]
            samp.phased = True
            continue

        if b in mom.alleles and b not in dad.alleles and a in dad.alleles:
            if fatherFirst:
                samp.alleles = [a, b]
            else:
                samp.alleles = [b, a]
            samp.phased = True
            continue
    
        if a in dad.alleles and a not in mom.alleles and b in mom.alleles:
            if fatherFirst:
                samp.alleles = [a, b]
            else:
                samp.alleles = [b, a]
            samp.phased = True
            continue

        if b in dad.alleles and b not in mom.alleles and a in mom.alleles:
            if fatherFirst:
                samp.alleles = [b, a]
            else:
                samp.alleles = [a, b]
            samp.phased = True
            continue

    return myVar

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # make a pedigree object from the PED file
    myPed = Pedigree(args.pedigree)

    # make Vcf object from the file
    myVcf = Vcf()

    # read and output the VCF header
    line = args.vcf.readline()
    while line and line[0] == '#':
        myVcf.add_header_info(line.rstrip())
        line = args.vcf.readline()
    myVcf.print_header()

    # parse each of the variants
    while line:
        myVar = Variant(myVcf, line.rstrip())
        s = myVar.sample_list[1]
        phasedVar = phase_by_transmission(myVar, myPed, args.fatherFirst, args.retain)
        print phasedVar.get_var_string()
        line = args.vcf.readline()

    # close the input file
    args.vcf.close()

# initialize the script
if __name__ == '__main__':
    sys.exit(main())
