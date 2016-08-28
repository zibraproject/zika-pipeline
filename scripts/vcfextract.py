#!/usr/bin/env python

import vcf
import sys
import subprocess
import csv
import os
from collections import defaultdict
from operator import attrgetter

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
                	depths[contig][int(pos)] = int(depth)
	return depths

for sample_tag in sys.argv[1:]:
	bamfn = "EM_079517_%s_marginalign.sorted.bam" % (sample_tag)
	vcffn = "np_EM_079517_%s.vcf" % (sample_tag)

	depths = collect_depths(bamfn)
	vcf_reader = vcf.Reader(filename=vcffn)
	for record in vcf_reader:
		print "%s\t%s\t%s\t%s" % (record.POS, record.QUAL, depths['EM_079517'][record.POS], record.INFO['TotalReads'][0])
