
import csv
import sys
from collections import defaultdict
import json

def get_runs_csv(fn):
	runs = defaultdict(defaultdict)
	reader = csv.reader(open(fn), delimiter='\t', dialect='excel-tab')
	lastrun = ''
	for row in reader:
		if row[11] != '':
			lastrun, notes = row[11].split('\n')

		if row[4] not in runs[lastrun]:
			runs[lastrun][row[4]] = list()
		runs[lastrun][row[4]].append(row[3])
	return runs

def get_runs():
    runs = defaultdict(defaultdict)
    runsamples = json.load(open('runsamples.json'))
    for sample in runsamples['data']:
        cols = sample['Run::run_name'].split("\n")
        run_name = cols[0]

        if not sample['sample_id'] in runs[run_name]:
            runs[run_name][sample['sample_id']] = list()

        runs[run_name][sample['sample_id']].append(sample['barcode_id'])
    return runs
