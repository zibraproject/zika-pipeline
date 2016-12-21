#!/bin/bash

ref=$1
dir=$2
scheme=$3
for barcode in NB01 NB02 NB03 NB04 NB05 NB06 NB07 NB08 NB09 NB10 NB11 NB12
do
  nanopolish extract -t 2d $dir/$barcode > $barcode.fastq
  variants.sh $ref $barcode $scheme $barcode
  margin_cons.py $ref $barcode.vcf $barcode.trimmed.sorted.bam a > $barcode.consensus.fasta
done

