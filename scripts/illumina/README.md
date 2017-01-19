#Analysis of Illumina Data

###Data Volume

Create a named data volume that mirrors local illumina data/ to data/ within container:

       docker create --name illumina-data -v /path/to/local/illumina/data/data:/illumina_data zibra/zibra

Modify align_reads/run_pipeline.sh and get_coverage/run_pipeline.sh to either execute snakemake on a single machine or using a scheduling system on a cluster.
Run align_reads/run_pipeline.sh to create alignments. To generate alignment statistics for merged bam files run get_coverage/run_pipeline.sh.
Output files are written to /build/illumina_analysis/ by default.