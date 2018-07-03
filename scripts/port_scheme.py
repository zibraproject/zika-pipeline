# update coordinates to match new reference

import sys
from Bio import SeqIO, Seq

primers = {}
newpos = {}

for ln in open(sys.argv[1], "r"):
	name, sequence = ln.rstrip().split("\t")
	primers[name] = sequence

recs = [rec for rec in SeqIO.parse(open(sys.argv[2]), "fasta")]
if len(recs) != 1:
	raise SystemExit

rec = recs[0]

for name, sequence in primers.items():
	pos = rec.seq.find(sequence)
	if pos != -1:
		newpos[name] = (pos, pos+len(sequence))
	else:
		rc = Seq.Seq(sequence).reverse_complement()
		pos = rec.seq.find(rc)
		if pos != -1:
			newpos[name] = (pos, pos+len(sequence))
			print ("Found %s revcom at %d" % (name, pos))
		else:
			print ("Not found %s" % (name,))


for ln in open(sys.argv[3], "r"):
	cols = ln.rstrip().split("\t")
	cols[0] = rec.id
	if cols[3] in newpos:
		cols[1] = str(newpos[cols[3]][0])
		cols[2] = str(newpos[cols[3]][1])
	print ("\t".join(cols))
