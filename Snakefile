shell.executable("/bin/bash")
shell.prefix("source $HOME/.bashrc; ")

IDS,  = glob_wildcards("runs/{id}.txt")
IDS2, = glob_wildcards("runs/{id}.txt")

rule all:
	input: expand("bam/{sample}.{sample2}.bam", sample=IDS, sample2=IDS2)

rule download:
	input: "runs/{id}.txt"
	output: 
		R1="fastq/{id}_1.fastq.gz",
		R2="fastq/{id}_2.fastq.gz"
	params:
		id="{id}"
	shell: "perl scripts/download.pl {params.id} && mv {params.id}*.fastq.gz fastq"


rule cutadapt:
	input:
		R1="fastq/{id}_1.fastq.gz",
		R2="fastq/{id}_2.fastq.gz"
	output:
		R1="trimmed/{id}_1.t.fastq.gz",
		R2="trimmed/{id}_2.t.fastq.gz"
	conda: "envs/cutadapt.yaml"
	shell: "cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o {output.R1} -p {output.R2} -O 5 --minimum-length=50 {input.R1} {input.R2}"


rule megahit:
	input:
		R1="trimmed/{id}_1.t.fastq.gz",
		R2="trimmed/{id}_2.t.fastq.gz"
	output: 
		di="megahit/{id}/",
		fa="megahit/{id}/final.contigs.fa"
	conda: "envs/megahit.yaml"
	threads: 8
	shell: "megahit --continue --k-list 27,47,67,87 --kmin-1pass -m 12e+10 --presets meta-large --min-contig-len 1000 -t {threads} -1 {input.R1} -2 {input.R2} -o {output.di}"


rule bwa_index:
	input:  "megahit/{id}/final.contigs.fa"
	output: "bwa_indices/{id}.fa.bwt"
	params:
		idx="bwa_indices/{id}.fa"
	conda: "envs/bwa.yaml"
	shell:
		'''
		bwa index -p {params.idx} {input}
		'''

rule bwa_mem:
	input:
		R1="trimmed/{id}_1.t.fastq.gz",
		R2="trimmed/{id}_2.t.fastq.gz",
		idx="bwa_indices/{id2}.fa.bwt"
	output: 
		bam="bam/{id}.{id2}.bam"
	params:
		idx="bwa_indices/{id}.fa"
	conda: "envs/bwa.yaml"
	shell: 
		'''
		bwa mem -t 4 {params.idx} {input.R1} {input.R2} | samtools sort -@4 -m 500M -o {output.bam} -
		samtools index {output.bam}

		samtools flagstat {output.bam} > {output.bam}.flagstat
		'''
	

