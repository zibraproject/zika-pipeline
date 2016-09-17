
import csv
import sys
from collections import defaultdict
#from runs import get_runs
import json

runs = defaultdict(defaultdict)

runsamples = json.load(open('runsamples.json'))
for sample in runsamples['data']:
    if not sample['sample_id'] in runs[sample['Run::run_name']]:
        runs[sample['Run::run_name']][sample['sample_id']] = list()

    runs[sample['Run::run_name']][sample['sample_id']].append(sample['barcode_id'])

for run, samples in runs.iteritems():
	print "cd %s" % (run,)
	#print "go_zika.sh ../newdata/%s/downloads/pass %s" % (run, sys.argv[2])
	for sample, barcodes in samples.iteritems():
		print "variants.sh %s %s %s" % (sample, sys.argv[1], " ".join(barcodes))
	print "cd .."
