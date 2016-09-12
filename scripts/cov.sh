#!/bin/bash

dir=$1

echo "| sample | reads | mapped | cov>=1 | cov>25 | link |"
echo "|--------|-------|--------|--------|--------|------|"
for barcode in "${@}"
do
    reads=`samtools view $barcode.sorted.bam| wc -l`
	mapped=`samtools view -F 4 $barcode.sorted.bam| wc -l`
    count1x=`samtools depth $barcode.trimmed.sorted.bam | awk '($3>0)' | wc -l`
    count25x=`samtools depth $barcode.trimmed.sorted.bam | awk '($3>=20)' | wc -l`
    echo "|" $barcode " | " $reads " | " $mapped " | " $count1x " | " $count25x " | "
	#<a href=\"http://s3.climb.ac.uk/nanopore/"$dir"_$barcode.tar\">FAST5 mapped</a> | "
done

