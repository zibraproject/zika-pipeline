#!/bin/bash -x
set -e

ref=$1
sample=$2
amplicons=$3

nanopolish extract --type 2d /data > $sample.fasta

# r9.4 models
cp -f /zibra/models/new_fast_template.5mers.model .
cp -f /zibra/models/new_fast_template.model .
cp -f /zibra/models/new_models.fofn .

bwa index $ref

bwa mem -xont2d $ref $sample.fasta | samtools view -bS - | samtools sort -o $sample.sorted.bam -
samtools index $sample.sorted.bam

align_trim.py --start --normalise 100 $amplicons --report $sample.alignreport.txt < $sample.sorted.bam 2> $sample.alignreport.er | samtools view -bS - | samtools sort -T $sample - -o $sample.trimmed.sorted.bam

align_trim.py --normalise 100 $amplicons --report $sample.alignreport.txt < $sample.sorted.bam 2> $sample.alignreport.er | samtools view -bS - | samtools sort -T $sample - -o $sample.primertrimmed.sorted.bam

samtools index $sample.trimmed.sorted.bam
samtools index $sample.primertrimmed.sorted.bam

#covplot.R $sample.alignreport.txt

nanopolish eventalign -vvvv -t 16 --models-fofn new_models.fofn --reads $sample.fasta -b $sample.trimmed.sorted.bam -g $ref --sam | samtools view -bS - | samtools sort -T $sample.tmp - -o $sample.np.sorted.bam
samtools index $sample.np.sorted.bam

nanopolish variants --progress -t 16 --reads $sample.fasta -o $sample.vcf -b $sample.trimmed.sorted.bam -e $sample.np.sorted.bam -g $ref -vv -w "`nanopolish_header.py $ref`" --snps --models-fofn new_models.fofn

nanopolish variants --progress -t 16 --reads $sample.fasta -o $sample.primertrimmed.vcf -b $sample.primertrimmed.sorted.bam -e $sample.np.sorted.bam -g $ref -vv -w "`nanopolish_header.py $ref`" --snps --models-fofn new_models.fofn

margin_cons.py $ref $sample.vcf $sample.trimmed.sorted.bam a > $sample.consensus.fasta
