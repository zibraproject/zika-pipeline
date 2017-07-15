#!/usr/bin/env python

#Written by Nick Loman (@pathogenomenick)
#Thanks to Aaron Quinlan for the argparse implementation from poretools.

import sys
import md5
import hashlib
import re
import argparse
import sqlite3
import version

def run_subtool(parser, args):
	if args.command == 'extract':
		import extract as submodule
	if args.command == 'basecaller':
		import basecaller as submodule
	if args.command == 'demultiplex':
		import demultiplex as submodule
	if args.command == 'minion':
		import minion as submodule

	# run the chosen submodule.
	submodule.run(parser, args)

class ArgumentParserWithDefaults(argparse.ArgumentParser):
	def __init__(self, *args, **kwargs):
		super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
		self.add_argument("-q", "--quiet", help="Do not output warnings to stderr",
						  action="store_true",
						  dest="quiet")

def main():
	parser = argparse.ArgumentParser(prog='zibra', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-v", "--version", help="Installed zibra version",
                        action="version",
                        version="%(prog)s " + str(version.__version__))
	subparsers = parser.add_subparsers(title='[sub-commands]', dest='command', parser_class=ArgumentParserWithDefaults)

	# extract
	parser_extract = subparsers.add_parser('extract',
                                          help='Create an empty poredb database')
	parser_extract.add_argument('directory', metavar='directory',
                             help='The name of the database.')
	parser_extract.add_argument('basecaller', metavar='basecaller',
                             help='The name of the basecaller')
	parser_extract.set_defaults(func=run_subtool)

	# callers
	parser_extract = subparsers.add_parser('basecaller', help='Display basecallers in files')
	parser_extract.add_argument('directory', metavar='directory', help='Directory of FAST5 files.')
	parser_extract.set_defaults(func=run_subtool)

	# demultiplex
	parser_demultiplex = subparsers.add_parser('demultiplex', help='Run demultiplex')
	parser_demultiplex.add_argument('fasta', metavar='fasta', help='Undemultiplexed FASTA file.')
	parser_demultiplex.add_argument('--threads', type=int, default=8, help='Number of threads')
	parser_demultiplex.add_argument('--prefix', help='Prefix for demultiplexed files')
	parser_demultiplex.set_defaults(func=run_subtool)

	# minion
	parser_minion = subparsers.add_parser('minion', help='Run demultiplex')
	parser_minion.add_argument('scheme', metavar='scheme', help='The name of the scheme.')
	parser_minion.add_argument('sample', metavar='sample', help='The name of the sample.')
	parser_minion.add_argument('--threads', type=int, default=8, help='Number of threads')
	parser_minion.add_argument('--scheme-directory', metavar='scheme_directory', default='/zibra/zika-pipeline/schemes', help='Default scheme directory')
	parser_minion.add_argument('--max-haplotypes', type=int, default=1000000, metavar='max_haplotypes', help='max-haplotypes value for nanopolish')
	parser_minion.set_defaults(func=run_subtool)

	# import
	"""parser_import = subparsers.add_parser('import',
	                                      help='Import files into a poredb database')
	parser_import.add_argument('db', metavar='DB',
	                           help='The poredb database.')
	parser_import.add_argument('fofn', metavar='FOFN',
	                           help='A file containing a list of file names.')
	parser_import.add_argument('--alternate-path', metavar='alternate_path')
	parser_import.set_defaults(func=run_subtool)
	"""

	args = parser.parse_args()

	if args.quiet:
		logger.setLevel(logging.ERROR)

	args.func(parser, args)

if __name__ == "__main__":
	main()
