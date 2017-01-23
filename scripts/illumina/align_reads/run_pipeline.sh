#!/bin/bash
#PBS -l walltime=120:00:00 -l nodes=1 -l mem=2gb -q new -o /gpfs/home/gkarthik/logs/snakelog.txt -j oe

cd ../db/
bwa index /zibra/zika-pipeline/scripts/illumina/db/zika_dc_2016.fa
cd ../align_reads/
# If running on laptop
snakemake
# If running on a cluster
# snakemake -j 50 --cluster-config cluster.json --cluster "qsub -l walltime={cluster.walltime} -l mem={cluster.mem} -l nodes={cluster.n} -q {cluster.queue} -o {cluster.logfile} -j {cluster.stdout}"
