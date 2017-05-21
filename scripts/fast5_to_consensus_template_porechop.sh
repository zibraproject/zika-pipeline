#!/bin/bash -x
set -e

schema=$1
sample=$2
directory=$3
flowcell=$4
barcode=$5

#nanopolish extract -r --type template ${directory} > ${sample}.fasta
#poretools fasta --type fwd --basecaller "ONT Albacore Sequencing Software=1.0.4" ${directory} > ${sample}.fasta
#python fasta.py ${directory} ${flowcell} > ${sample}.pretrim.fasta 2> ${sample}.skipped.txt
#porechop --untrimmed -i ${sample}.pretrim.fasta -b porechop-${sample} --barcode_threshold 75 --threads 16 --check_reads 1000 --barcode_diff 2 --require_two_barcodes

cp porechop-${directory}/NB${barcode}.fasta.gz ${sample}.fasta.gz
gunzip -f ${sample}.fasta.gz

ref=/zibra/zika-pipeline/schemes/${schema}/V1/${schema}.reference.fasta
bed=/zibra/zika-pipeline/schemes/${schema}/V1/${schema}.scheme.bed

# Takes a $sample.fasta file full of nanopolish reads, ie
# nanopolish extract --type 2d /data > $sample.fasta
# Files are written to working directory

# 2) copy the r9.4 model files into current directory
#cp -f /zibra/models/new_fast_template.5mers.model .
#cp -f /zibra/models/new_fast_template.model .
#cp -f /zibra/models/new_models.fofn .

# 3) index the ref & align with bwa
bwa index $ref
bwa mem -xont2d $ref $sample.fasta | samtools view -bS - | samtools sort -o $sample.sorted.bam -
samtools index $sample.sorted.bam

# 4) trim the alignments to the primer start sites and normalise the coverage to save time
align_trim.py --start --normalise 100 $bed --report $sample.alignreport.txt < $sample.sorted.bam 2> $sample.alignreport.er | samtools view -bS - | samtools sort -T $sample - -o $sample.trimmed.sorted.bam
align_trim.py --normalise 100 $bed --report $sample.alignreport.txt < $sample.sorted.bam 2> $sample.alignreport.er | samtools view -bS - | samtools sort -T $sample - -o $sample.primertrimmed.sorted.bam
samtools index $sample.trimmed.sorted.bam
samtools index $sample.primertrimmed.sorted.bam

#covplot.R $sample.alignreport.txt

# 6) do variant calling using the raw signal alignment
nanopolish variants --progress -t 16 --reads $sample.fasta -o $sample.vcf -b $sample.trimmed.sorted.bam -g $ref -vv -w "`nanopolish_header.py $ref`" --snps --ploidy 1
nanopolish variants --progress -t 16 --reads $sample.fasta -o $sample.primertrimmed.vcf -b $sample.primertrimmed.sorted.bam -g $ref -vv -w "`nanopolish_header.py $ref`" --snps --ploidy 1

# 7) variant frequency plot
vcfextract.py ${sample} > ${sample}.variants.tab

# 8) filter the variants and produce a consensus
margin_cons.py $ref $sample.vcf $sample.trimmed.sorted.bam a > $sample.consensus.fasta
