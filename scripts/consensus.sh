#!/bin/bash -x

ref_prefix=$1
sample=$2
poretools_dir=$3
sample_tag=$4
second_batch=$5

#python vcftagdepth.py np_"$ref_prefix"_"$sample_tag".vcf "$ref_prefix"_"$sample_tag"_marginAlign.sorted.bam > "$sample_tag"_"$ref_prefix"_marginAlign.tagged.vcf
#python vcffilter_pp.py "$sample_tag"_"$ref_prefix"_marginAlign.tagged.vcf > "$sample_tag"_"$ref_prefix"_marginAlign_0.3.tagged.vcf
#python vcftagnp.py "$sample_tag"_"$ref_prefix"_marginAlign.tagged.vcf np_"$ref_prefix"_"$sample_tag".vcf > "$sample_tag"_"$ref_prefix"_np.tagged.vcf
vcftagprimersites.py all np_"$ref_prefix"_"$sample_tag".vcf > "$sample_tag"_"$ref_prefix"_np_primer.tagged.vcf
vcffilter.py "$sample_tag"_"$ref_prefix"_np_primer.tagged.vcf > "$sample_tag"_"$ref_prefix"_np_primer.filtered075_30.vcf
vcffilterqual.py "$sample_tag"_"$ref_prefix"_np_primer.tagged.vcf > "$sample_tag"_"$ref_prefix"_np_primer.filtered_qual200.vcf
#cp "$sample_tag"_"$ref_prefix"_np_primer.filtered.vcf ../vcfs
#python margin_cons.py "$ref_prefix".fasta ../vcfs/"$sample_tag"_"$ref_prefix"_np_primer.filtered.vcf ../bams/"$ref_prefix"_"$sample_tag"_marginAlign.sorted.bam  2>"$ref_prefix"_"$sample_tag".stderr | python tagfasta.py ../samples.txt "$sample" > "$ref_prefix"_"$sample_tag".consensus.fasta
#bedtools genomecov -d -ibam ../bams/"$ref_prefix"_"$sample_tag"_marginAlign.sorted.bam > "$ref_prefix"_"$sample_tag"_marginAlign.sorted.bam.cov
#Rscript covplot.R "$ref_prefix"_"$sample_tag".bedtools_19rx_bwa.txt "$ref_prefix"_"$sample_tag".cov.pdf
#Rscript covplot.R "$ref_prefix"_"$sample_tag".bedtools_19rx_bwa.txt "$ref_prefix"_"$sample_tag".cov.png
#coverageBed -abam ../bams/"$ref_prefix"_"$sample_tag"_marginAlign.sorted.bam -b ../11_rx_v1.bed -d | groupBy -g 1,2,3,4 -c 6 -o mean > "$ref_prefix"_"$sample_tag".bedtools.txt
