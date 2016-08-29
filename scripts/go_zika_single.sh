#!/bin/bash

prefix=$1
dir=$2

poretools fastq --type 2D $dir > $prefix.fastq
variants.sh $prefix
margin_cons.py refs/Zika_FP.fasta $barcode.vcf $barcode.trimmed.sorted.bam a > $barcode.consensus.fasta
