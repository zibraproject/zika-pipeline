#!/bin/bash -x

ref_prefix=$1
sample=$2
poretools_dir=$3
sample_tag=$4
second_batch=$5
read_type=$6

tmp=`tempfile`
poretools fastq --type 2D ../data/$poretools_dir/pass > $tmp
if [ "$read_type" == "hq" ];
then
	poretools fastq --type 2D --high-quality ../data/$poretools_dir/fail >> $tmp
elif [ "$read_type" == "pass" ]
then
	echo "nothing"
fi

if [ ! "$second_batch" == "na" ]
then
	poretools fastq --type 2D ../data/$second_batch/pass >> $tmp
	if [ "$read_type" == "hq" ];
	then
        	poretools fastq --type 2D --high-quality ../data/$second_batch/fail >> $tmp
	elif [ "$read_type" == "pass" ]
	then
        	echo "nothing"
	fi
	poretools fastq --type 2D ../data/$second_batch/pass >> $tmp
fi
mv $tmp "$sample_tag".fastq

rm -rf jobTree_align_"$ref_prefix"_"$sample_tag"

#marginAlign "$sample_tag".fastq ../refs/"$ref_prefix".fasta "$ref_prefix"_"$sample_tag".sam --em --outputModel "$sample_tag".output.hmm --jobTree jobTree_align_"$ref_prefix"_"$sample_tag" --maxThreads 40
marginAlign "$sample_tag".fastq ../refs/"$ref_prefix".fasta "$ref_prefix"_"$sample_tag".sam --jobTree jobTree_align_"$ref_prefix"_"$sample_tag" --maxThreads 40 --inputModel input.hmm
samtools view -bS "$ref_prefix"_"$sample_tag".sam | samtools sort - "$ref_prefix"_"$sample_tag"_marginalign.sorted
samtools index "$ref_prefix"_"$sample_tag"_marginalign.sorted.bam
