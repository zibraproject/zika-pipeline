
import csv
import sys
from collections import defaultdict

runs = defaultdict(defaultdict)
reader = csv.reader(open(sys.argv[1]), delimiter='\t', dialect='excel-tab')
lastrun = ''
for row in reader:
	if row[9] != '':
		lastrun = row[9]

	if row[4] not in runs[lastrun]:
		runs[lastrun][row[4]] = list()
	runs[lastrun][row[4]].append(row[3])

for run, samplelookup in runs.iteritems():
	print "# %s" % ( run , )

	for s, barcodes in samplelookup.iteritems():
		print "variants.sh %s ../zika-pipeline/metadata/v2.amplicons.bed %s" % (s, " ".join(barcodes))
		print "margin_cons.py ../refs/Zika_FP.fasta %s.vcf %s.trimmed.sorted.bam a > %s.consensus.fasta" % (s, s, s)
