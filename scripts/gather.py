import sys
from Bio import SeqIO
import tempfile
import os
import glob
import shutil

# extract with constraints:
#   -- only one group ever
#   -- only one flowcell ID ever
#   -- always unique read ID

def run(parser, args):
	d = '%s/workspace/pass' % (args.directory)

	all_fastq_outfn = "%s_all.fastq" % (os.path.split(args.directory)[-1])
	all_fastq_outfh = open(all_fastq_outfn, "w")

	for root, dirs, files in os.walk(d):
		paths = os.path.split(root)
		barcode_directory = paths[-1]

		fastq = [root+'/'+f for f in files if f.endswith('.fastq')]
		if len(fastq):
			fastq_outfn = "%s_%s.fastq" % (os.path.split(args.directory)[-1], barcode_directory)
			outfh = open(fastq_outfn, "w")
			print >>sys.stderr, "Processing %s files in %s" % (len(fastq), barcode_directory)

			dups = set()
			uniq = 0
			total = 0	
			for f in fastq:
				for rec in SeqIO.parse(open(f), "fastq"):
					total += 1
					if rec.id not in dups:
						SeqIO.write([rec], outfh, "fastq")
						SeqIO.write([rec], all_fastq_outfh, "fastq")

						dups.add(rec.id)
						uniq += 1

			outfh.close()

			print "%s\t%s\t%s" % (fastq_outfn, total, uniq)

	all_fastq_outfh.close()

