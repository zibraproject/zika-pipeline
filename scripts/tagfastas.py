#!/usr/bin/env python

import sqlite3
import re
import sys
from Bio import SeqIO

con = sqlite3.connect(sys.argv[1])
con.row_factory = sqlite3.Row
cur = con.cursor()

def lookup_sample(sample):
	cur.execute("select * from samples, runs where runs.Batch = ? and runs.sample_fk = samples.rowid", (sample,))
	row = cur.fetchone()
	return row

for rec in SeqIO.parse(sys.stdin, "fasta"):
	m = re.search(r'EM_079517_(.*)_hq', rec.id)
	sample = m.group(1)

	metadata = lookup_sample(sample)

	if metadata['prefecture'] == 'Kambia':
		country = 'SLE'
	else:
		country = 'GUI' 
	loc = '-'.join([metadata['prefecture'], metadata['sousprefecture'], metadata['village'].encode('ascii', 'ignore')])
	rec.id = '|'.join(('EBOV', metadata['LabID'], 'MinION', country, loc, metadata['date_sample_taken'])).replace(' ', '_')
	rec.description = ""
	SeqIO.write([rec], sys.stdout, "fasta")

