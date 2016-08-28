#!/usr/bin/env python

import pysam
import sys
from copy import copy
from vcftagprimersites import read_bed_file

def trim(s, start_pos, end):
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

bed = read_bed_file('all')

def find_primer(pos, direction):
	# {'Amplicon_size': '1874', 'end': 7651, '#Region': 'region_4', 'start': 7633, 'Coords': '7633', "Sequence_(5-3')": 'GCTGGCCCGAAATATGGT', 'Primer_ID': '16_R'}
	from operator import itemgetter

	closest = min([(abs(p['start'] - pos), p['start'] - pos, p) for p in bed if p['direction'] == direction], key=itemgetter(0))
	return closest

infile = pysam.AlignmentFile("-", "rb")
outfile = pysam.AlignmentFile("-", "wh", template=infile)
for s in infile:
	cigar = copy(s.cigartuples)

	if len(sys.argv) > 1:
		if not s.query_name.startswith(sys.argv[1]):
			continue

	## logic - if alignment start site is _before_ but within X bases of
	## a primer site, trim it off

	if s.is_unmapped:
		continue

	p1 = find_primer(s.reference_start, 'F')
	p2 = find_primer(s.reference_end, 'R')

	print >>sys.stderr, "%s\t%s\t%s_%s\t%s\t%s\t%s\t%s" % (s.reference_start, s.reference_end, p1[2]['Primer_ID'], p2[2]['Primer_ID'], p1[2]['Primer_ID'], abs(p1[1]), p2[2]['Primer_ID'], abs(p2[1]))

	#amp_len = p1['end'] - p2['end']

	## if the alignment starts before the end of the primer, trim to that position
	#print >>sys.stderr, s.reference_start, p1[2]['end']
	#print >>sys.stderr, s.reference_end, p2[2]['start']

	"""
	if p1[0] < 50:
		if s.reference_start < p1[2]['start'] - 1:
	#		trim(s, p1[2]['end']-1, 0)
			trim(s, p1[2]['start']-1, 0)
	else:
		trim(s, s.reference_start + 20, 0)

	if p2[0] < 50:
		if s.reference_end > p2[2]['end'] - 1:
	#		trim(s, p2[2]['start']-1, 1)
			trim(s, p2[2]['end']-1, 1)
	else:
		trim(s, s.reference_end - 20, 1)
	"""

	try:
		trim(s, s.reference_start + 20, 0)
		trim(s, s.reference_end - 20, 1)

		outfile.write(s)
	except Exception:
		pass

	outfile.write(s)
