# Docker notes

## Build

From within the `zika-pipeline` directory. Build `zibra` image:

    docker build -t zibra/zibra:latest .

Push image to Docker Hub:

    docker push zibra/zibra:latest

## Run

### Data volume

Create a named data volume that mirrors local `data/` to `data/` within container:

    docker create --name zibra-data -v /Volumes/Meristem/data:/data zibra/zibra    

This is to get data into the Docker container. Note that the path to local directory has to be an absolute path.

Create a named data volume for a single sample:

    docker create --name zibra-data-lb01-nb01 -v /Volumes/Meristem/data/usvi-library1-2016-12-10/basecalled_reads/pass_demultiplex/NB01:/data zibra/zibra

### Build volume

Create a named data volume that mirrors local `build/` to `build/` within container:

    docker create --name zibra-build -v /Users/trvrb/Documents/src/zika-pipeline/build:/build zibra/zibra

This is to get data out of the Docker container. Note that the path to local directory has to be an absolute path.

### Start

Enter docker image:

    docker run -t -i --volumes-from zibra-data --volumes-from zibra-build zibra/zibra /bin/bash

Run single sample script within image:

    ./scripts/go_single_sample_r94.sh refs/KJ776791.2.fasta NB03 metadata/v2_500.amplicons.ver2.bed
