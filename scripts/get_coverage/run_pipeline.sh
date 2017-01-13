#!/bin/bash
#PBS -l walltime=120:00:00 -l nodes=1 -q new -o /gpfs/home/gkarthik/logs/new.txt -j oe

cd /gpfs/home/gkarthik/jobs/bam_operations/get_coverage/
snakemake -j 50 --cluster-config cluster.json --cluster "qsub -l walltime={cluster.walltime} -l mem={cluster.mem} -l nodes={cluster.n} -q {cluster.queue} -o {cluster.logfile} -j {cluster.stdout}"
