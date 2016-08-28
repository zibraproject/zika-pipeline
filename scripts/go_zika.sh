#!/bin/bash

dir=$1
for barcode in NB01 NB02 NB03 NB04 NB05 NB06 NB07 NB08 NB09 NB10 NB11 NB12
do
  poretools fastq --type 2D $dir/$barcode > $barcode.fastq
  variants.sh $barcode
  margin_cons.py ../refs/Zika_FP.fasta $barcode.vcf $barcode.trimmed.sorted.bam a > $barcode.consensus.fasta
done

