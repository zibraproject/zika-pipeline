#!/bin/bash -x

ref_prefix=$1
sample=$2
poretools_dir=$3
sample_tag=$4
second_batch=$5
read_type=$6

align_trim.py < "$ref_prefix"_"$sample_tag".sam >/dev/null 2>"$ref_prefix"_"$sample_tag".alignments.txt
