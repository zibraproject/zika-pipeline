from poretools.Fast5File import Fast5FileSet
import sys
import tempfile
import os
import shutil

# extract with constraints:
#   -- only one group ever
#   -- only one flowcell ID ever
#   -- always unique read ID

def run(parser, args):
	tmpdir = tempfile.mkdtemp(dir='.')

	cmd = ("porechop --verbosity 2 --untrimmed -i \"%s\" -b %s --barcode_threshold 80 --threads %s --check_reads 10000 --barcode_diff 5 --require_two_barcodes > %s.demultiplexreport.txt" % (args.fasta, tmpdir, args.threads, args.fasta))
	print >>sys.stderr, cmd
	os.system(cmd)

	a, b = os.path.split(args.fasta)
	prefix, ext = os.path.splitext(b)

	for fn in os.listdir(tmpdir):
		newfn = "%s-%s" % (prefix, os.path.basename(fn))
		shutil.move(tmpdir + '/' + fn, newfn)

		os.system("gunzip -f %s" % (newfn,))

	if not args.no_remove_directory:
		os.rmdir(tmpdir)

