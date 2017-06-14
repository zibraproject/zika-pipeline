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

	cmd = ("porechop --untrimmed -i \"%s\" -b %s --barcode_threshold 75 --threads %s --check_reads 1000 --barcode_diff 2 --require_two_barcodes" % (args.fasta, tmpdir, args.threads))
	print >>sys.stderr, cmd
	os.system(cmd)

	a, b = os.path.split(args.fasta)
	prefix, ext = os.path.splitext(b)

	for fn in os.listdir(tmpdir):
		newfn = "%s-%s" % (prefix, os.path.basename(fn))
		shutil.move(tmpdir + '/' + fn, newfn)

		os.system("gunzip -f %s" % (newfn,))

	os.rmdir(tmpdir)

