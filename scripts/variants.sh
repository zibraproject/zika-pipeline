#!/bin/bash -x
set -e

ref=$1
sample=$2
amplicons=$3
barcode1=$4
barcode2=$5

seqtk seq -A "$barcode1".fastq >$sample.fasta

if [ -n "$barcode2" ]; then
   seqtk seq -A "$barcode2".fastq >>$sample.fasta
fi

bwa index $ref

bwa mem -xont2d $ref $sample.fasta | samtools view -bS - | samtools sort -o $sample.sorted.bam -
samtools index $sample.sorted.bam

align_trim.py --start --normalise 100 $amplicons --report $sample.alignreport.txt < $sample.sorted.bam 2> $sample.alignreport.er | samtools view -bS - | samtools sort -T $sample - -o $sample.trimmed.sorted.bam

align_trim.py --normalise 100 $amplicons --report $sample.alignreport.txt < $sample.sorted.bam 2> $sample.alignreport.er | samtools view -bS - | samtools sort -T $sample - -o $sample.primertrimmed.sorted.bam

samtools index $sample.trimmed.sorted.bam
samtools index $sample.primertrimmed.sorted.bam

covplot.R $sample.alignreport.txt

nanopolish eventalign -vvvv -t 16 --models-fofn ../models/nanopolish_models.fofn --reads $sample.fasta -b $sample.trimmed.sorted.bam -g $ref --sam | samtools view -bS - | samtools sort -T $sample.tmp - -o $sample.np.sorted.bam
samtools index $sample.np.sorted.bam

nanopolish variants --progress -t 16 --reads $sample.fasta -o $sample.vcf -b $sample.trimmed.sorted.bam -e $sample.np.sorted.bam -g $ref -vv -w "`nanopolish_header.py $ref`" --snps --models-fofn ../models/nanopolish_models.fofn

nanopolish variants --progress -t 16 --reads $sample.fasta -o $sample.primertrimmed.vcf -b $sample.primertrimmed.sorted.bam -e $sample.np.sorted.bam -g $ref -vv -w "`nanopolish_header.py $ref`" --snps --models-fofn ../models/nanopolish_models.fofn


