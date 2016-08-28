#!/bin/bash -x
set -e

ref_prefix=$1
sample=$2
poretools_dir=$3
sample_tag=$4
second_batch=$5
read_type=$6

tmp=`tempfile`
poretools fastq --type 2D ../data/$poretools_dir/pass > $tmp
poretools fasta --type 2D ../data/$poretools_dir/pass > "$tmp".fasta
if [ "$read_type" == "hq" ];
then
	poretools fastq --type 2D --high-quality ../data/$poretools_dir/fail >> $tmp
	poretools fasta --type 2D --high-quality ../data/$poretools_dir/fail >> "$tmp".fasta
elif [ "$read_type" == "pass" ]
then
	echo "nothing"
fi

if [ ! "$second_batch" == "na" ]
then
	poretools fastq --type 2D ../data/$second_batch/pass >> $tmp
	poretools fasta --type 2D ../data/$second_batch/pass >> "$tmp".fasta
	if [ "$read_type" == "hq" ];
	then
        	poretools fastq --type 2D --high-quality ../data/$second_batch/fail >> $tmp
        	poretools fasta --type 2D --high-quality ../data/$second_batch/fail >> "$tmp".fasta
	elif [ "$read_type" == "pass" ]
	then
        	echo "nothing"
	fi
fi
#head -20000 $tmp > $tmp.20000
#mv $tmp.20000 "$sample_tag".fasta
sed 's/_Basecall_2D_2d/_Basecall_2D_000_2d/' $tmp > "$sample_tag".fastq
sed 's/_Basecall_2D_2d/_Basecall_2D_000_2d/' "$tmp".fasta > "$sample_tag".fasta

if [ ! -e ../refs/"$ref_prefix".fasta.amb ]
then
    bwa index ../refs/"$ref_prefix".fasta
fi
#bwa mem -t 1 -x ont2d ../refs/"$ref_prefix".fasta "$sample_tag".fasta  | samtools view -bS - | samtools sort - "$ref_prefix"_"$sample_tag"_bwa.sorted
#coverageBed -abam "$ref_prefix"_"$sample_tag"_bwa.sorted.bam -b ../metadata/11_rx_v3.bed -d | groupBy -g 1,2,3,4 -c 6 -o mean > "$ref_prefix"_"$sample_tag".bedtools_11rx_bwa.txt
#coverageBed -abam "$ref_prefix"_"$sample_tag"_bwa.sorted.bam -b ../metadata/19_rx.bed -d | groupBy -g 1,2,3,4 -c 6 -o mean > "$ref_prefix"_"$sample_tag".bedtools_19rx_bwa.txt
#samtools index "$ref_prefix"_"$sample_tag"_bwa.sorted.bam

rm -rf jobTree_align_"$ref_prefix"_"$sample_tag"
marginAlign "$sample_tag".fastq ../refs/"$ref_prefix".fasta "$ref_prefix"_"$sample_tag".sam --jobTree jobTree_align_"$ref_prefix"_"$sample_tag" --inputModel ../models/input_10_04.hmm
align_trim.py < "$ref_prefix"_"$sample_tag".sam | samtools view -bS - | samtools sort - "$ref_prefix"_"$sample_tag"_marginalign.sorted
samtools index "$ref_prefix"_"$sample_tag"_marginalign.sorted.bam

rm -f "$sample_tag".fasta.fast5.fofn
nanopolish eventalign --reads "$sample_tag".fasta -b "$ref_prefix"_"$sample_tag"_marginalign.sorted.bam -g ../refs/"$ref_prefix".fasta --sam | samtools view -bS - | samtools sort - "$ref_prefix"_"$sample_tag"_np.sorted
samtools index "$ref_prefix"_"$sample_tag"_np.sorted.bam
nanopolish variants --models-fofn offset_models.fofn --progress -t 1 --reads "$sample_tag".fasta -o np_"$ref_prefix"_"$sample_tag".vcf -b "$ref_prefix"_"$sample_tag"_marginalign.sorted.bam -e "$ref_prefix"_"$sample_tag"_np.sorted.bam -g ../refs/"$ref_prefix".fasta -vv -w "EM_079517:0-20000" --snp

