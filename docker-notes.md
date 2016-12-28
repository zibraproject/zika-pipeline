# Docker notes

## Build

From within the `zika-pipeline` directory. Build `zibra` image:

    docker build -t zibra/zibra:latest .

Push image to Docker Hub:

    docker push zibra/zibra:latest

## Run

### Data

Create a named data volume that mirrors local `zika-pipeline/data/` to `/data/` within container:

    docker create --name zibra-data -v /Users/trvrb/Documents/src/zika-pipeline/data:/data zibra/zibra    

_Note that the path to local directory has to be an absolute path._

Create a named data volume for a single sample:

    docker create --name zibra-data-lb01-nb03 -v /Users/trvrb/Documents/src/zika-pipeline/data/libraries/usvi-library1-2016-12-10/basecalled_reads/pass_demultiplex/NB03:/data zibra/zibra  

### Start

Enter docker image:

    docker run -t -i --volumes-from zibra-data zibra/zibra /bin/bash

Run single sample script within image:

    ./scripts/go_single_sample_r94.sh refs/KJ776791.2.fasta NB03 metadata/v2_500.amplicons.ver2.bed
