#!/bin/bash -x
set -e

ref_prefix=$1
sample=$2
poretools_dir=$3
sample_tag=$4
second_batch=$5
read_type=$6

#tmp=`tempfile`
#poretools fastq --type 2D $poretools_dir > $tmp
#poretools fasta --type 2D $poretools_dir > "$tmp".fasta
#mv $tmp "$sample_tag".fastq
#mv "$tmp".fasta "$sample_tag".fasta

if [ ! -e ../refs/"$ref_prefix".fasta.amb ]
then
    bwa index ../refs/"$ref_prefix".fasta
fi

rm -rf jobTree_align_"$ref_prefix"_"$sample_tag"
marginAlign "$sample_tag".fastq ../refs/"$ref_prefix".fasta "$ref_prefix"_"$sample_tag".sam --jobTree jobTree_align_"$ref_prefix"_"$sample_tag" --inputModel ../models/input_10_04.hmm
align_trim.py < "$ref_prefix"_"$sample_tag".sam | samtools view -bS - | samtools sort -T "$ref_prefix"_"$sample_tag" - -o "$ref_prefix"_"$sample_tag"_marginalign.sorted.bam
cat "$ref_prefix"_"$sample_tag".sam | samtools view -bS - | samtools sort -T "$ref_prefix"_"$sample_tag" - -o "$ref_prefix"_"$sample_tag"_marginalign.sorted.bam
samtools index "$ref_prefix"_"$sample_tag"_marginalign.sorted.bam

rm -f "$sample_tag".fasta.fast5.fofn
nanopolish eventalign --reads "$sample_tag".fasta -b "$ref_prefix"_"$sample_tag"_marginalign.sorted.bam -g ../refs/"$ref_prefix".fasta --sam | samtools view -bS - | samtools sort -T "$ref_prefix"_"$sample_tag"_np - -o "$ref_prefix"_"$sample_tag"_np.sorted.bam
samtools index "$ref_prefix"_"$sample_tag"_np.sorted.bam
nanopolish variants --progress -t 1 --reads "$sample_tag".fasta -o np_"$ref_prefix"_"$sample_tag".vcf -b "$ref_prefix"_"$sample_tag"_marginalign.sorted.bam -e "$ref_prefix"_"$sample_tag"_np.sorted.bam -g ../refs/"$ref_prefix".fasta -vv -w "gi|992324757|gb|KU707826.1|:0-11000" --snps

