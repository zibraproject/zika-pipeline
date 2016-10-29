#!/usr/bin/env python

import vcf
import sys
import subprocess
import csv
import os
from collections import defaultdict
from operator import attrgetter
from runs import get_runs

def read_vcf(fn):
    vcfinfo = {}
    vcf_reader = vcf.Reader(open(fn, 'r'))
    for record in vcf_reader:
        vcfinfo[record.POS] = record
    return vcfinfo

def collect_depths(bamfile):
    if not os.path.exists(bamfile):
        raise SystemExit("bamfile %s doesn't exist" % (bamfile,))

    print >>sys.stderr, bamfile

    p = subprocess.Popen(['samtools', 'depth', bamfile],
                             stdout=subprocess.PIPE)
    out, err = p.communicate()
    depths = defaultdict(dict)
    for ln in out.split("\n"):
            if ln:
                    contig, pos, depth = ln.split("\t")
                    depths[int(pos)] = int(depth)
    return depths

positions = set()

runs = get_runs()
for run, samples in runs.iteritems():
    for sample_tag in samples.keys():
        for vcfset in ['', '.primertrimmed']:
            vcffn = "%s/%s%s.vcf" % (run, sample_tag, vcfset)
            if not os.path.exists(vcffn):
                continue

            print >>sys.stderr, vcffn

            vcf_reader = vcf.Reader(filename=vcffn)
            for record in vcf_reader:
                if len(record.ALT[0]) == 1 and len(record.REF) == 1:
                    positions.add(record.POS)


print "pos\tset\trun\tsample\tvartype\tdepth\tsupportfraction\tbasecalledfrequency\tsuppportingreads"

for run, samples in runs.iteritems():
    for sample_tag in samples.keys():
        for vcfset in ['', '.primertrimmed']:
            vcffn = "%s/%s%s.vcf" % (run, sample_tag, vcfset)
            if not os.path.exists(vcffn):
                continue

            vcffile = read_vcf(vcffn)
            bamfn = "%s/%s.primertrimmed.sorted.bam" % (run, sample_tag)
            depths = collect_depths(bamfn)

            #1-based pyvcf
            for pos in positions:
                if pos-1 in depths:
                    depth = depths[pos-1]
                else:
                    depth = 0

                if pos in vcffile:
                    info = vcffile[pos].INFO
                    print "%s\t%s\t%s\t%s\tvariant\t%s\t%s\t%s\t%s" % (pos, vcfset, run, sample_tag, depth, info['SupportFraction'][0], info['BaseCalledFrequency'][0], info['SupportingReads'][0])
                else:
                    print "%s\t%s\t%s\t%s\tinvariant\t%s\t0\t0\t0" % (pos, vcfset, run, sample_tag, depth)

