#!/bin/bash

sample=$1
fq1=$1
fq2=$2
scheme=$3

snakemake -p --verbose --config scheme="${scheme}" fq1="${fq1}" fq2="${fq2}" sample="${sample}"  --snakefile /zibra/zika-pipeline/scripts/illumina/align_reads/Snakefile all

