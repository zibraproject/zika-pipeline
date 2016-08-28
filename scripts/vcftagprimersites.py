#!/usr/bin/env python

import vcf
import sys
import subprocess
import csv
from collections import defaultdict

sets = {
	'all'  : '../metadata/all_primers.txt',
	'19rx' : '../metadata/19_rx.txt',
	'11rx' : '../metadata/11_rx_v3.txt'
}

def read_bed_file(primerset):
	sets = {
		'all'  : '../metadata/all_primers.txt',
	    '19rx' : '../metadata/19_rx.txt',
	    '11rx' : '../metadata/11_rx_v3.txt'
	}

	bedfile = []
	with open(sets[primerset]) as csvfile:
		reader = csv.DictReader(csvfile, dialect='excel-tab')
		for row in reader:
			row['direction'] = row['Primer_ID'].split("_")[1]
			row['end'] = int(row['Coords']) + len(row["Sequence_(5-3')"])
			row['start'] = int(row['Coords'])
			bedfile.append(row)
	return bedfile

def overlaps(coords, pos):
	for v in coords:
		if pos >= v['start'] and pos <= v['end']:
			return v
	return False

if __name__ == "__main__":
	if sys.argv[1] not in sets:
		print "Invalid set"
		raise SystemExit

	bedfile = read_bed_file(sys.argv[1])

	vcf_reader = vcf.Reader(filename=sys.argv[2])
	vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
	for record in vcf_reader:
		v = overlaps(bedfile, record.POS)
		if v:
			record.INFO['PRIMER'] = v["Sequence_(5-3')"]

#	PP = list(record.INFO)
#	record.INFO = {}
#	record.INFO['PP'] = PP
#	record.INFO['DEPTH'] = depths[record.CHROM][record.POS]

		vcf_writer.write_record(record)
