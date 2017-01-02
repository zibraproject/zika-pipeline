# Zika bioinformatics pipeline

Run scripts via [Docker](https://www.docker.com/). Notes on building and pushing Docker image are [here](docker-notes.md).

## Run

### Data volume

Create a named data volume that mirrors local MinION `data/` to `data/` within container:

    docker create --name zibra-data -v /Volumes/Meristem/data:/data zibra/zibra

This is to get MinION data into the Docker container. Note that the path to local directory has to be an absolute path. Change `/Volumes/Meristem/data` to wherever local data is stored. [Notes on data schema are here](data-schema.md).

### Samples volume

Create a named data volume that mirrors local sample metadata `samples/` to `samples/` within container:

    docker create --name zibra-samples -v /Volumes/Meristem/samples:/samples zibra/zibra

This is to get sample metadata into the Docker container. Note that the path to local directory has to be an absolute path. Change `/Volumes/Meristem/samples` to wherever local data is stored. [Notes on metadata schema are here](data-schema.md)

### Build volume

Create a named data volume that mirrors local `build/` to `build/` within container:

    docker create --name zibra-build -v /Volumes/Meristem/build:/build zibra/zibra

This is to get data out of the Docker container. Note that the path to local directory has to be an absolute path. Change `/Volumes/Meristem/build` to wherever local results are stored.

### Start

Enter docker image:

    docker run -t -i --volumes-from zibra-data --volumes-from zibra-samples --volumes-from zibra-build zibra/zibra /bin/bash

Run script:

    python scripts/pipeline.py
