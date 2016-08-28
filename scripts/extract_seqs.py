#!/usr/bin/env python
from Bio import SeqIO
import sys

for rec in SeqIO.parse(sys.stdin, "fasta"):
	if rec.id in sys.argv[1:]:
		SeqIO.write([rec], sys.stdout, "fasta")
