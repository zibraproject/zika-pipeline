# Docker notes

## Build

From within the `zika-pipeline` directory. Build `zibra` image:

    docker build -t zibra/zibra:latest .

Push image to Docker Hub:

    docker push zibra/zibra:latest

## Clean

Remove all containers:

    docker rm `docker ps --no-trunc -aq`
