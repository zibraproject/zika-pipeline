#!/usr/bin/env python

import pysam
import sys
from copy import copy
from vcftagprimersites import read_bed_file
from collections import defaultdict

def trim(cigar, s, start_pos, end):
	if not end:
		pos = s.pos
	else:
		pos = s.reference_end

	eaten = 0
	while 1:
		## chomp stuff off until we reach pos
		if end:
			flag, length = cigar.pop()
		else:
			flag, length = cigar.pop(0)

		#print >>sys.stderr,  "Chomped a %s, %s" % (flag, length)
		if flag == 0:
			## match
			#to_trim -= length
			eaten += length
			if not end:
				pos += length
			else:
				pos -= length
		if flag == 1:
			## insertion to the ref
			#to_trim -= length
			eaten += length
		if flag == 2:
			## deletion to the ref
			#eaten += length
			if not end:
				pos += length
			else:
				pos -= length
			pass
		if flag == 4:
			eaten += length
		if not end and pos >= start_pos and flag == 0:
			break
		if end and pos <= start_pos and flag == 0:
			break

	#print >>sys.stderr, "pos:%s %s" % (pos, start_pos)

	extra = abs(pos - start_pos)
	#print >> sys.stderr, "extra %s" % (extra)
	if extra:
		if flag == 0:
			#print >>sys.stderr,  "Inserted a %s, %s" % (0, extra)

			if end:
				cigar.append((0, extra))
			else:
				cigar.insert(0, (0, extra))
			eaten -= extra

	if not end:
			s.pos = pos - extra

	#print >>sys.stderr,  "New pos: %s" % (s.pos)

	if end:
		cigar.append((4, eaten))
	else:
		cigar.insert(0, (4, eaten))
	oldcigarstring = s.cigarstring
	s.cigartuples = cigar

	#print >>sys.stderr,  s.query_name, oldcigarstring[0:50], s.cigarstring[0:50]

def find_primer(bed, pos, direction):
	# {'Amplicon_size': '1874', 'end': 7651, '#Region': 'region_4', 'start': 7633, 'Coords': '7633', "Sequence_(5-3')": 'GCTGGCCCGAAATATGGT', 'Primer_ID': '16_R'}
	from operator import itemgetter

	closest = min([(abs(p['start'] - pos), p['start'] - pos, p) for p in bed if p['direction'] == direction], key=itemgetter(0))
	return closest


def go(args):
    bed = read_bed_file(args.bedfile)

    counter = defaultdict(int)

    infile = pysam.AlignmentFile("-", "rb")
    outfile = pysam.AlignmentFile("-", "wh", template=infile)
    for s in infile:
            cigar = copy(s.cigartuples)

            ## logic - if alignment start site is _before_ but within X bases of
            ## a primer site, trim it off

            if s.is_unmapped:
                    continue

            if s.is_supplementary:
                    continue

            p1 = find_primer(bed, s.reference_start, '+')
            p2 = find_primer(bed, s.reference_end, '-')

            print >>sys.stderr, "%s\t%s\t%s\t%s_%s\t%s\t%s\t%s\t%s\t%s\t%s" % (s.query_name, s.reference_start, s.reference_end, p1[2]['Primer_ID'], p2[2]['Primer_ID'], p1[2]['Primer_ID'], abs(p1[1]), p2[2]['Primer_ID'], abs(p2[1]), s.is_secondary, s.is_supplementary)

            ## if the alignment starts before the end of the primer, trim to that position

            try:
                if s.reference_start < p1[2]['end']:
                    trim(cigar, s, p1[2]['end'], 0)
                else:
                    continue

                if s.reference_end > p2[2]['end']:
                    trim(cigar, s, p2[2]['end'], 1)
                else:
                    continue
            except Exception, e:
                print >>sys.stderr, "problem %s" % (e,)
                pass

            if args.normalise:
                pair = "%s-%s-%d" % (p1[2]['Primer_ID'], p2[2]['Primer_ID'], s.is_reverse)
                counter[pair] += 1

                if counter[pair] > args.normalise:
                    continue

            ## if the alignment starts before the end of the primer, trim to that position

#	try:
#		trim(s, s.reference_start + 40, 0)
#		trim(s, s.reference_end - 40, 1)
#
#		outfile.write(s)
#	except Exception:
#		pass

            outfile.write(s)


import argparse

parser = argparse.ArgumentParser(description='Trim alignments from an amplicon scheme.')
parser.add_argument('bedfile', help='BED file containing the amplicon scheme')
parser.add_argument('--normalise', type=int, help='Subsample to n coverage')

args = parser.parse_args()
go(args)

