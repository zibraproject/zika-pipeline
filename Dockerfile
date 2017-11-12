# Zibra pipeline Docker container
FROM ubuntu:16.04
MAINTAINER Nick Loman <n.j.loman@bham.ac.uk>

# base software
RUN apt-get update
RUN apt-get install -y git build-essential wget zlib1g-dev vim libncurses5-dev

# python and python dependencies
RUN apt-get install -y python python-pip
RUN pip install pysam pyvcf biopython clint

# Add Python3 for snakemake and plot_coverage.py for illumina run
RUN apt-get install -y python3 python3-pip
RUN pip3 install snakemake pandas

# create working directory
RUN mkdir /zibra
WORKDIR /zibra/

# BWA
RUN git clone --recursive https://github.com/lh3/bwa && cd bwa && make -j8

# samtools
RUN apt-get install -y libbz2-dev liblzma-dev && wget https://github.com/samtools/samtools/releases/download/1.4/samtools-1.4.tar.bz2 && tar xvjf samtools-1.4.tar.bz2 && cd samtools-1.4 && make

# Smith-Waterman library
RUN git clone https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library.git && cd Complete-Striped-Smith-Waterman-Library/src && make

# Poretools
RUN pip install git+https://github.com/arq5x/poretools.git@basecaller-choice

# porechop branch
RUN pip3 install git+https://github.com/rrwick/Porechop.git

# version - cache bust
ADD HISTORY /zibra/HISTORY

# nanopolish
RUN git clone --recursive https://github.com/jts/nanopolish/ && cd nanopolish && git checkout tags/v0.8.4 && make -j4

# zibra pipeline
RUN git clone https://github.com/zibraproject/zika-pipeline
WORKDIR /zibra/zika-pipeline/

# environmental variables
ENV PATH $PATH:/zibra/nanopolish:/zibra/bwa:/zibra/samtools-1.4:/zibra/zika-pipeline/scripts
ENV PYTHONPATH /zibra/Complete-Striped-Smith-Waterman-Library/src
ENV LD_LIBRARY_PATH /zibra/Complete-Striped-Smith-Waterman-Library/src
