# GoogleMAGs
Workflow for MAG construction

## Snakemake arguments
$ snakemake  -p --verbose --keep-remote  -j [number of available cores] --kubernetes -s Snakefile --default-remote-provider GS  --default-remote-prefix [GCS Bucket] --use-conda --container-image '[Path to Docker image if different than default]'
Examlpe

$ snakemake  -p --verbose --keep-remote  -j 400 --kubernetes -s Snakefilev11 --default-remote-provider GS  --default-remote-prefix hn-snakemake --use-conda --container-image 'gcr.io/tagareby/snakemake'


## Building a custom Docker Image on the 
You can start from the sample Dockerfile in this repo which just addes procrps which includes the linux tool free required for megahit -m flag when value less than 1

- Auth docker using gcloud
$ gcloud auth configure-docker
- Build the local DOckerfile
$ docker build -t [user]/snakemake .
- Tag the local image
$ docker tag [user]/snakemake gcr.io/[project_name]/snakemake
- Push the image to the Google Cloud Container Registery
$ docker push gcr.io/[project]/snakemake
