#!/bin/bash -x
prefix=$1
newprefix=$2
threshold=$3


#make tree for each cluster
cat "$prefix"_aligned.fasta | shorten-headers.py > "$prefix"_aligned_short.fasta
rm ${newprefix}*cluster*
cp "$prefix"_aligned_short.fasta "$newprefix"_aligned_short.fasta
clustertree.R "$prefix"_aligned_trim_short_raxml.nwk "$threshold" | tee "$newprefix"_clusters.txt | split-clusters.py "$newprefix"_aligned_short.fasta

for i in `ls ${newprefix}*cluster[0-9]`; 
do
 n=`echo "$i" | awk '{print substr($0,length-1,2)}'`
 if [ "$n" != 'r0' ]; then
   cat ../refs/earliest.fasta | shorten-headers.py >> ${i}
 fi

 rm RAxML_*"${i}"_raxml*;
 raxmlHPC-PTHREADS-SSE3 -T 4 -f a -x 12345 -p 12345 -# 100 -m GTRGAMMA -s ${i} -n "${i}"_raxml;
 cat ${i} | get-alignment.py 2>${i}.positions > ${i}.alignment; 

 #make pdf
  xvfb-run pdf_tree.py --positions ${i}.positions --alignment ${i}.alignment RAxML_bipartitions.${i}_raxml "$prefix"_microreact_private.csv ${i}.pdf
 sleep 5;
  xvfb-run pdf_tree.py --positions ${i}.positions --alignment ${i}.alignment RAxML_bipartitions.${i}_raxml "$prefix"_microreact_private.csv ${i}.png
 sleep 5;

done



