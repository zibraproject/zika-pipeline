#!/usr/bin/env python

from Bio import SeqIO
import sys

SeqIO.convert(sys.argv[1], "fastq", sys.argv[2], "fasta")
