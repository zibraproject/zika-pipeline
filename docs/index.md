# The ZiBRA Pipeline

During the Zika project we have developed a fully integrated sequencing pipeline which is described in our paper at Nature Protocols.

This page describes in more detail the bioinformatics pipeline from that paper.

## Installing Docker

For Linux:

https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-16-04

For Windows:

Ideally Docker for Windows:

https://docs.docker.com/docker-for-windows/

Some users may need to use Docker Toolbox:

https://docs.docker.com/toolbox/toolbox_install_windows/

For Mac:

https://docs.docker.com/docker-for-mac/

## Quick start

### MinION pipeline

Install the Docker container:

```
  docker pull zibra/zibra:latest
  docker run -t -i zibra/zibra:latest /bin/bash
```

Download and run on the Zika test sample:

```
  mkdir whoref
  cd whoref
  wget https://s3.climb.ac.uk/nanopore/Zika_Control_Material_R9.4_2D.tar
  tar xvf Zika_Control_Material_R9.4_2D.tar
  fast5_to_consensus.sh ZikaAsian WHO 20161118_Zika/downloads/pass/NB08
```

A number of files will be produced of interest:



### Analysing diverse viral lineages

We designed this protocol in order to monitor evolutionary changes over short time-scales for example in outbreaks and epidemics that have been through a recent bottleneck. In practice it means this protocol has been mainly tested on virus genomes with <1% divergence from a reference, such as during the recent Ebola and Zika epidemics. The pipeline relies on detecting differences from the reference genome, that in turn depends on the nanopolish software. Nanopolish works by testing all possible combinations of variants to determine which combination is best explained by the signal-level data produced by the nanopore device. In theory it is possible to work with more divergent viral lineages, but computational complexity rapidly increases as the sequence diverges (2^n, where n is the number of variants to test in a window).

By default the Zibra pipeline will test a maximum of 1,000,000 potential haplotypes in a window (typically 100 bp), which roughly equates to 2^20 mutations. This can be increased using the --max-haplotypes parameter to zibra.py, but at the expense of running time.

### Local offline basecalling

Currently we are recommending Oxford Nanopore’s Albacore for local base calling. If basecalling offline on Windows or Mac, we recommend this is done on the native operating system rather than in the Docker container, for reasons of speed.

We also do not recommend using Albacore demultiplexing at present, as this is rather lenient as it requires only a single barcode copy.

To basecall an R9.4 1D dataset on Windows with eight threads use, if the reads are stored in ``C:\data\reads\run`` and you want the basecals to save to ``C:\Users\nick\data\run`` then run::

``read_fast5_basecaller.py --input C:\data\reads\run --worker_threads 8 -c r94_450bps_linear.cfg -s C:\Users\nick\data\run -r -o fast5``

### Demultiplexing

Demultiplexing (splitting reads by their barcodes) is currently rather convoluted in the current Zibra pipeline. Initially, reads are mapped against a reference genome and only mapped reads are considered further. These reads are then trimmed by considering their alignment coordinates and relating that to the known start and end coordinates of each amplicon (as determined by the specified Primal Scheme BED file), with 40 base flanking sequences retained.

![Trimming](trimming.png)

Next, Porechop is run: the presence of double barcodes are enforced, with a minimum barcode identity of 80%. Because nanopolish cannot read trimmed files, the results of Porechop are used to reconstitute the original FASTA files but this time sorted into bins, and then the process is repeated to call consensus sequences.

The reason for this complicated process is to reduce any cross-barcode contamination caused through the presence of in silico chimeras, e.g. when the start or end of a read is erroneously detected. This is important for accurate determination of contamination in the negative control.

### Mounting a local directory

You will probably want to mount a local directory where your reads are kept.

On Windows, say my data is in ``C:\Users\nick\data`` you would run:

docker run -t -v //C/Users/nick/data:/data -i zibra/zibra:latest /bin/bash

On Mac:

docker run -t -v /Users/nick/data:/data -i zibra/zibra:latest /bin/bash

Then use Porechop (in the container) to demultiplex the reads:

### Illumina pipeline- Quickstart

   docker pull zibra/zibra:latest
   docker run -t -i zibra/zibra:latest /bin/bashwget --no-check-certificate 
   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/007/SRR5122847/SRR5122847_1.fastq.gz
   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/007/SRR5122847/SRR5122847_2.fastq.gz
   illumina_pipeline.sh SRR5122847 SRR5122847_1.fastq.gz SRR5122847_2.fastq.gz ZikaAsian

## Credits

The ZiBRA Pipeline was developed with contributions from:

  - Nick Loman (MinION pipeline)
  - Jared Simpson (nanopolish SNP calling)
  - Matt Loose (nanopore demultiplexing script)
  - Karthik Gangavarapu (Illumina pipeline)
  - Nate Grubaugh (Illumina pipeline)
  - Kristian Andersen (Illumina pipeline)
  - Trevor Bedford (help with Docker and useful fixes)

It relies on a whole heap of open source software, thank you to all contributors.



