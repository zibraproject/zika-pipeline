# Docker notes

## Build

From within the `zika-pipeline` directory. Build `zibra` image:

    docker build -t zibra/zibra:latest .

Push image to Docker Hub:

    docker push zibra/zibra:latest

## Run

### Data

Create a named data volume that mirrors local `zika-pipeline/data/` to `data/` within container:

    docker create --name zibra-data -v /Users/trvrb/Documents/src/zika-pipeline/data:/data zibra/zibra    

_Note that the path to local directory has to be an absolute path._

### Start

Run docker image:

    docker run -t -i --volumes-from zibra-data zibra/zibra /bin/bash
