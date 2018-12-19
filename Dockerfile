#Sample custom snakemake container based on the official image
FROM quay.io/snakemake/snakemake
RUN apt-get update && apt-get install -y procps && apt-get clean -y
