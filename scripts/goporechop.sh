#!/bin/bash -x
set -e

library=$1
flowcell=$2

python fasta.py /data/"$library"_albacore/workspace $flowcell > $library.fasta
porechop --untrimmed -i "$library".fasta -b porechop-"$library" --barcode_threshold 75 --threads 16 --check_reads 1000 --barcode_diff 2 --require_two_barcodes
