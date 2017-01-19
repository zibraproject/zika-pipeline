# Zibra pipeline Docker container
FROM ubuntu:16.04
MAINTAINER Nick Loman <n.j.loman@bham.ac.uk>

# base software
RUN apt-get update
RUN apt-get install -y git build-essential wget zlib1g-dev vim libncurses5-dev

# python and python dependencies
RUN apt-get install -y python python-pip
RUN pip install pysam pyvcf biopython

# Add Python3 for snakemake and plot_coverage.py for illumina run
RUN apt-get install -y python3 python3-pip
RUN pip3 install snakemake
RUN pip3 install pandas

# create working directory
RUN mkdir /zibra
WORKDIR /zibra/

# nanopolish
RUN git clone --recursive https://github.com/jts/nanopolish.git && cd nanopolish && make -j8

# BWA
RUN git clone --recursive https://github.com/lh3/bwa && cd bwa && make -j8

# samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && tar xvjf samtools-1.3.1.tar.bz2 && cd samtools-1.3.1 && make

# R9.4 models for nanopolish
RUN mkdir models && cd models && wget http://s3.climb.ac.uk/nanopore/nanopolish_r94models.tar && tar xvf nanopolish_r94models.tar

# Smith-Waterman library
RUN git clone https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library.git && cd Complete-Striped-Smith-Waterman-Library/src && make

# zibra pipeline
RUN git clone https://github.com/zibraproject/zika-pipeline
WORKDIR /zibra/zika-pipeline/

# update
ADD http://www.timeapi.org/utc/now /tmp/bustcache
RUN git pull

# environmental variables
ENV PATH $PATH:/zibra/nanopolish:/zibra/bwa:/zibra/samtools-1.3.1:/zibra/zika-pipeline/scripts
ENV PYTHONPATH /zibra/Complete-Striped-Smith-Waterman-Library/src
ENV LD_LIBRARY_PATH /zibra/Complete-Striped-Smith-Waterman-Library/src