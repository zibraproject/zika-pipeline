#!/bin/bash

dir=$1
sample=$2

poretools fasta --type 2D $dir > $sample.fasta
demultiplex.py --barcodes /zibra/zika-pipeline/barcodes/barcodes.fasta $sample.fasta
