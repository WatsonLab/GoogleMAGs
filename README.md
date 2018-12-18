# GoogleMAGs
Workflow for MAG construction

##Snakemake arguments
$ snakemake  -p --verbose --keep-remote  -j [number of available cores] --kubernetes -s Snakefile --default-remote-provider GS  --default-remote-prefix [GCS Bucket] --use-conda --container-image '[Path to Docker image if different than default]'
Examlpe

$ snakemake  -p --verbose --keep-remote  -j 400 --kubernetes -s Snakefilev11 --default-remote-provider GS  --default-remote-prefix hn-snakemake --use-conda --container-image 'gcr.io/tagareby/snakemake'
