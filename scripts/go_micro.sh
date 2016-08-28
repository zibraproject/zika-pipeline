#!/bin/bash -x
consensus=$1
prefix=$2
threads=$3

#make large tree

#if [ ! -a "$consensus" ]
#then
#  echo "$consensus does not exist!"
#  exit
#fi

cat $consensus ../refs/earliest.fasta ../refs/related_sle_isolates.fasta | sed 's/ /_/' > "$prefix"_aligned.fasta
trim.py "$prefix"_aligned.fasta 50 18500 > "$prefix"_aligned_trim.fasta
# FastTree -nt "$prefix"_aligned_trim.fasta > "$prefix"_aligned_trim_fasttree.nwk

raxmlHPC-PTHREADS-SSE3 -T $threads -f a -x 12345 -p 12345 -# 100 -m GTRGAMMA -s "$prefix"_aligned_trim.fasta -n "$prefix"_raxml

# root it

root_and_deheader.py RAxML_bipartitions."$prefix"_raxml "EBOV|EM_079422|EMLab|GUI|Macenta--|2014-03-27" > "$prefix"_aligned_trim_short_raxml.nwk

#make microreact csv
generate_csv_ig.py ../metadata/prefectures.txt "$prefix"_aligned_trim.fasta public > "$prefix"_microreact_public.csv
generate_csv_ig.py ../metadata/prefectures.txt "$prefix"_aligned_trim.fasta private > "$prefix"_microreact_private.csv

