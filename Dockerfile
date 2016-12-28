# Zibra pipeline Docker container
FROM ubuntu:16.04
MAINTAINER Nick Loman <n.j.loman@bham.ac.uk>
RUN apt-get update && apt-get install -y \
                      git \
					  build-essential \
					  wget \
					  zlib1g-dev \
					  vim \
					  libncurses5-dev \
					  python \
					  python-pip
# Python dependencies
RUN pip install pysam pyvcf biopython
# nanopolish
RUN mkdir /zibra && cd /zibra/ && git clone --recursive https://github.com/jts/nanopolish.git && cd nanopolish && make -j8
# BWA
RUN cd /zibra/ && git clone --recursive https://github.com/lh3/bwa && cd bwa && make -j8
# Zibra pipeline scripts
RUN cd /zibra/ && git clone https://github.com/zibraproject/zika-pipeline
# Samtools
RUN cd /zibra && wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && tar xvjf samtools-1.3.1.tar.bz2 && cd samtools-1.3.1 && make
# R9.4 models for nanopolish
RUN cd /zibra && mkdir models && cd models && wget http://s3.climb.ac.uk/nanopore/nanopolish_r94models.tar && tar xvf nanopolish_r94models.tar
# Environmental variables
ENV PATH $PATH:/zibra/nanopolish:/zibra/bwa:/zibra/samtools-1.3.1:/zibra/zika-pipeline/scripts
