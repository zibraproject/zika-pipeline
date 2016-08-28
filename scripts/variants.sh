#!/bin/bash -x
set -e

sample=$1

seqtk seq -A $sample.fastq > $sample.fasta

bwa mem -xont2d ../refs/Zika_FP.fasta $sample.fastq | samtools view -bS - | samtools sort -o $sample.sorted.bam -
samtools index $sample.sorted.bam

align_trim.py < $sample.sorted.bam | samtools view -bS - | samtools sort -T $sample - -o $sample.trimmed.sorted.bam
samtools index $sample.trimmed.sorted.bam

nanopolish eventalign -t 16 --models-fofn ../models/reftrained.modelset.fofn --reads $sample.fasta -b $sample.trimmed.sorted.bam -g ../refs/Zika_FP.fasta --sam | samtools view -bS - | samtools sort -T $sample.tmp - -o $sample.np.sorted.bam
samtools index $sample.np.sorted.bam

nanopolish variants --progress -t 16 --reads $sample.fasta -o $sample.vcf -b $sample.trimmed.sorted.bam -e $sample.np.sorted.bam -g ../refs/Zika_FP.fasta -vv -w "gi|631250742|gb|KJ776791.1|:0-11000" --snps --models-fofn ../models/reftrained.modelset.fofn


