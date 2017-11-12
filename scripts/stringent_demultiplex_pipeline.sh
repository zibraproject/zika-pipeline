#!/bin/bash

tag=$1
scheme=$2
zibra.py extract /data/data/"$tag"_albacore/workspace/ > "$tag".fasta
zibra.py minion --skip-nanopolish --normalise 0 $scheme "$tag"
python /zibra/zika-pipeline/scripts/align_trim_fasta.py "$tag".trimmed.sorted.bam /zibra/zika-pipeline/schemes/"$scheme"/V1/"$scheme".scheme.bed > "$tag".trimmed.fasta 2> "$tag".trimmed.report.txt
zibra.py demultiplex --threads 96 "$tag".trimmed.fasta
for bc in 01 02 03 04 05 06 07 08 09 10 11 12; do python /zibra/zika-pipeline/scripts/reconstitute.py "$tag".fasta "$tag".trimmed-BC"$bc".fasta > "$tag"-BC"$bc"_originals.fasta; done

