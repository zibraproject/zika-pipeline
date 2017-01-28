#!/bin/bash
#PBS -l walltime=120:00:00 -l nodes=1 -q new -o /gpfs/home/gkarthik/logs/new.txt -j oe

snakemake
# If executing on a cluster
# snakemake -j 50 --cluster-config cluster.json --cluster "qsub -l walltime={cluster.walltime} -l mem={cluster.mem} -l nodes={cluster.n} -q {cluster.queue} -o {cluster.logfile} -j {cluster.stdout}"
