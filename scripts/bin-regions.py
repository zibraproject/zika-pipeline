#!/usr/bin/env python

## read a bam file

## look at start & end position
## bin to closest F and R pair
## print out

## tag       start primer   end primer      count 
## sample1   1_F            2_R             X

import sys
import sqlite3

con = sqlite3.connect(sys.argv[1])
con.row_factory = sqlite3.Row
cur = con.cursor()

def get_runs(dataset):
    if dataset == 'all':
        cur.execute("select samples.ROWID, * from runs, samples where runs.sample_fk = samples.ROWID")
    else:
        cur.execute("select samples.ROWID, * from runs, samples where runs.sample_fk = samples.ROWID and runs.dataset = ?", (dataset,))
    return cur.fetchall()

runs = get_runs(sys.argv[2])
for r in runs:
	fn = "EM_079517_%s_hq.alignments.txt" % (r['Batch'])
	for ln in open(fn, "r"):
		print r['LabID'] + "\t" + ln,	

