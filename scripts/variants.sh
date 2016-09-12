#!/bin/bash -x
set -e

sample=$1
amplicons=$2
barcode1=$3
barcode2=$4

seqtk seq -A "$barcode1".fastq >$sample.fasta

if [ -n "$barcode2" ]; then
   seqtk seq -A "$barcode2".fastq >>$sample.fasta
fi

bwa index ../refs/Zika_FP.fasta

bwa mem -xont2d ../refs/Zika_FP.fasta $sample.fasta | samtools view -bS - | samtools sort -o $sample.sorted.bam -
samtools index $sample.sorted.bam

align_trim.py --normalise 100 $amplicons --report $sample.alignreport.txt < $sample.sorted.bam 2> $sample.alignreport.er | samtools view -bS - | samtools sort -T $sample - -o $sample.trimmed.sorted.bam
samtools index $sample.trimmed.sorted.bam

covplot.R $sample.alignreport.txt

nanopolish eventalign -t 16 --models-fofn ../models/nanopolish_models.fofn --reads $sample.fasta -b $sample.trimmed.sorted.bam -g ../refs/Zika_FP.fasta --sam | samtools view -bS - | samtools sort -T $sample.tmp - -o $sample.np.sorted.bam
samtools index $sample.np.sorted.bam

nanopolish variants --progress -t 16 --reads $sample.fasta -o $sample.vcf -b $sample.trimmed.sorted.bam -e $sample.np.sorted.bam -g ../refs/Zika_FP.fasta -vv -w "gi|631250742|gb|KJ776791.1|:0-11000" --snps --models-fofn ../models/nanopolish_models.fofn


