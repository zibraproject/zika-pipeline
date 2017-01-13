#!/bin/bash
#PBS -l walltime=120:00:00 -l nodes=1 -l mem=2gb -q new -o /gpfs/home/gkarthik/logs/snakelog.txt -j oe

cd /gpfs/home/gkarthik/jobs/zika-pipeline/
snakemake -j 50 --cluster-config cluster.json --cluster "qsub -l walltime={cluster.walltime} -l mem={cluster.mem} -l nodes={cluster.n} -q {cluster.queue} -o {cluster.logfile} -j {cluster.stdout}"
