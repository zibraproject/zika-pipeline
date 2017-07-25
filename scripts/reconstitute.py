#!/usr/bin/env python

#Written by Nick Loman
#Part of the ZiBRA pipeline (zibraproject.org)


import argparse
import sys
import glob
from Bio import SeqIO
import gzip

def go(args):
	ids = [rec.id for rec in SeqIO.parse(gzip.open(args.trimmed), "fasta")]
	SeqIO.write((seq for seq in SeqIO.parse(args.fasta, "fasta") if seq.id in ids), sys.stdout, "fasta")

parser = argparse.ArgumentParser(description='Trim alignments from an amplicon scheme.')
parser.add_argument('fasta', help='Original FASTA file')
parser.add_argument('trimmed', help='Trimmed FASTA GZ file')

args = parser.parse_args()
go(args)

