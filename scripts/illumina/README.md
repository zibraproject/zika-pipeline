#Analysis of Illumina Data

###Data Volume

Create a named data volume that mirrors the local .fastq.gz files generated to data/ within container:

       docker create --name illumina-data -v /path/to/local/illumina/data/data:/illumina_data zibra/zibra

Modify align_reads/run_pipeline.sh and get_coverage/run_pipeline.sh to execute snakemake on a single machine or using a scheduling system on a cluster. By default the script executes the pipeline on the local machine.

Run ./run_pipeline.sh in align_reads/ to create aligned bam files. Output files are written to a folder named using the timestamp in /build/illumina_analysis/ by default. To generate alignment statistics for merged bam files run ./run_pipeline.sh in the get_coverage/ directory. Modify the src parameter in get_coverage/Snakefile to point to the new output folder. The statistics file is written to a _pileup/ direcoty in _aligned_bams/ by default. 

Example directory structure of output folder,

```
.
|-- _aligned_bams
|   |-- Sample1.aligned.sorted.bam
|   |-- Sample1.trimmed.aligned.sorted.bam
|   |-- Sample2.aligned.sorted.bam
|   |-- Sample2.trimmed.aligned.sorted.bam
|   `-- _pileup
|       |-- Sample1.aligned.sorted.tsv
|       |-- Sample1.trimmed.aligned.sorted.tsv
|       |-- Sample2.aligned.sorted.tsv
|       |-- Sample2.trimmed.aligned.sorted.tsv
|       `-- statistics.md
|-- _reads
|   |-- Sample1_R1.fastq
|   |-- Sample1_R2.fastq
|   |-- Sample2_R1.fastq
|   `-- Sample2_R2.fastq
`-- _reports
    |-- Sample1.alignreport.txt
    `-- Sample1.alignreport.txt
```
