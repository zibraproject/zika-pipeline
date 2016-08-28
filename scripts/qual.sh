#!/bin/bash -x

ref_prefix=$1
sample=$2
poretools_dir=$3
sample_tag=$4
second_batch=$5
read_type=$6

#bwa mem -t 1 -x ont2d ../refs/"$ref_prefix".fasta "$sample_tag".fasta  | samtools view -bS - | samtools sort - "$ref_prefix"_"$sample_tag"_bwa.sorted
#expand-cigar.py --fasta ../refs/"$ref_prefix".fasta --bam "$ref_prefix"_"$sample_tag"_bwa.sorted.bam | count-errors.py - > "$ref_prefix"_"$sample_tag"_bwa.idystats.txt
expand-cigar.py --fasta ../refs/"$ref_prefix".fasta --bam "$ref_prefix"_"$sample_tag"_marginalign.sorted.bam | count-errors.py - > "$ref_prefix"_"$sample_tag"_marginalign.idystats.txt
