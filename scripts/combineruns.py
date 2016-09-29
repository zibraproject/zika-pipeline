
import csv
import sys
from collections import defaultdict
from runs import get_runs

runs = get_runs()

for run, samples in runs.iteritems():
	print "cd %s" % (run,)
	print "go_zika.sh ../newdata/%s/downloads/pass %s" % (run, sys.argv[1])
	for sample, barcodes in samples.iteritems():
		print "variants.sh %s %s %s" % (sample, sys.argv[1], " ".join(barcodes))
	print "cd .."
