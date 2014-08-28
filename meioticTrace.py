#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2014-01-30 09:11 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
meiotic trace.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Trace meiotic lineage of 3 generation pedigree")
    #parser.add_argument('-r', '--retain', action='store_true', help='retain original phasing (default: overwrite)')
    #parser.add_argument('-f', '--fatherFirst', action='store_true', help='outputs phasing as father|mother (default is mother|father)')
    parser.add_argument('-s', '--sample', type=str, required=True, help='sample to trace meiotic lineage. Must be phased and have family info in pedigree')
    parser.add_argument('-p', '--pedigree', type=argparse.FileType('r'), required=True, help='PED format file (http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped)') 
    parser.add_argument('vcf', nargs='?', type=argparse.FileType('r'), default=None, help='phased vcf file to read. Phasing must be maternal|paternal. If \'-\' or absent then defaults to stdin.')

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

    # returns whether variant has only two alleles: 0 and 1
    def is_diallelic(self):
        return len(self.alt) == 1

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

    def is_het(self):
        return len(self.alleles) == len(set(self.alleles))

    def is_hom_ref(self):
        return self.alleles == ['0','0']

    def is_hom_alt(self):
        return len(set(self.alleles)) == 1 and '0' not in self.alleles

    def is_hom(self):
        return len(set(self.alleles)) == 1

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

# check for mendelian transmission of alleles
def is_mendelian(myVar, pedigree, samp_id):
    samp_ped = pedigree.get_sample(samp_id)
    mom_id = samp_ped.mat
    dad_id = samp_ped.pat

    # set samp and parents
    samp = myVar.sample_info(samp_id)
    mom = myVar.sample_info(mom_id)
    dad = myVar.sample_info(dad_id)

    if len(samp.alleles) == 2 and '.' not in samp.alleles:
        a = samp.alleles[0]
        b = samp.alleles[1]
    else:
        return False

    # if a not in mom.alleles, then b must have come from the mom, and a must have come from the dad
    if a not in mom.alleles:
        if a not in dad.alleles or b not in mom.alleles:
            return False
    if b not in mom.alleles:
        if b not in dad.alleles or a not in mom.alleles:
            return False
    # same for dad
    if a not in dad.alleles:
        if a not in mom.alleles or b not in dad.alleles:
            return False
    if b not in dad.alleles:
        if b not in mom.alleles or a not in dad.alleles:
            return False

    #print "mendelian error: %s: %s, %s: %s, %s: %s" % (samp_id, samp.get_gt_string(), mom_id, mom.get_gt_string(), dad_id, dad.get_gt_string())
    
    # if you got to this point, then you're mendelian
    return True

# primary function
def meiotic_trace(myVar, pedigree, samp_id):
    # the output list [maternal allele origin, paternal allele origin]
    # the first item is the origin of the sample's maternal allele, and the
    # second item is the origin of the sample's paternal allele.
    # each item is a 2 digit binary number.
        # 00: allele came from mother's mother
        # 01: allele came from mother's father
        # 10: allele came from father's mother
        # 11: allele came from father's father
        # -1: unknown
    trace = [-1] * 2

    # bail if more than 1 alt allele
    if not myVar.is_diallelic():
        return trace

    #print samp_id
    samp_ped = pedigree.get_sample(samp_id)
    mom_id = samp_ped.mat
    dad_id = samp_ped.pat

    # set samp
    if samp_id in myVar.sample_list:
        samp = myVar.sample_info(samp_id)

    # if the sample is phased, store a and b as the alleles in the sample
    if samp.is_phased() and len(samp.alleles) == 2 and '.' not in samp.alleles:
        a = samp.alleles[0]
        b = samp.alleles[1]
    else:
        return trace

    # check that the mom and dad are in the VCF
    if mom_id in myVar.sample_list and dad_id in myVar.sample_list:
        mom = myVar.sample_info(mom_id)
        dad = myVar.sample_info(dad_id)
    else:
        return trace
    
    # skip if parent has nocall
    if not mom.is_phased() or not dad.is_phased() or '.' in mom.alleles + dad.alleles:
        return trace

    # check for mendelian transmission
    if not is_mendelian(myVar, pedigree, samp_id):
        return trace

    # only look at heterozygous sites in the sample where one of the parents in het
    if samp.is_hom() or (mom.is_hom() and dad.is_hom()):
        return trace

    # output is gonna be a list of length 2.
    # the first item is the origin of the sample's maternal allele, and the
    # second item is the origin of the sample's paternal allele.
    # each item is a 2 digit binary number.
        # 00: allele came from mother's mother
        # 01: allele came from mother's father
        # 10: allele came from father's mother
        # 11: allele came from father's father
        # -1: unknown

    # only informative when at least one of the parents is het
    if mom.is_het():
        if a == mom.alleles[0]:
            trace[0] = 0
        elif a == mom.alleles[1]:
            trace[0] = 1
    if dad.is_het():
        if b == dad.alleles[0]:
            trace[1] = 10
        elif b == dad.alleles[1]:
            trace[1] = 11

    return trace

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
    if line[0] != '#':
        print 'VCF header not found'
        exit(1)
    while line and line[0] == '#':
        myVcf.add_header_info(line.rstrip())
        line = args.vcf.readline()
    #myVcf.print_header()

    # parse each of the variants
    while line:
        myVar = Variant(myVcf, line.rstrip())
        s = myVar.sample_list[1]
        #phasedVar = phase_by_transmission(myVar, myPed, args.fatherFirst, args.retain)
        #print phasedVar.get_var_string()
        trace = meiotic_trace(myVar, myPed, args.sample)
        # if trace is uninformative [-1,-1] then skip
        if sum(trace) > -2:
            print '\t'.join(map(str,[myVar.chrom, myVar.pos, myVar.ref, ','.join(myVar.alt), myVar.qual] + trace ))
            #print '\t'.join(map(str,trace)) + '\t' + myVar.get_var_string()
        line = args.vcf.readline()

    # close the input file
    args.vcf.close()

# initialize the script
if __name__ == '__main__':
    sys.exit(main())
