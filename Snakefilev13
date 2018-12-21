#Adding support to GCS
from pathlib import Path
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()

GS_INPUT = "hn-snakemake/big"
GS_INPUT = "mw_test_metagenomics"
GS_PREFIX = GS_INPUT

shell.executable("/bin/bash")
shell.prefix("source $HOME/.bashrc; ")

IDS,  = glob_wildcards("runs/{id}.txt")
IDS2, = glob_wildcards("runs/{id}.txt")

rule all:
	input: GS.remote(expand(GS_PREFIX + "/coverage/{sample}.{sample2}.txt", sample=IDS, sample2=IDS2))

#rule download has been integrated into the cutadpat to reduce network transfer

rule cutadapt:
	input: "runs/{id}.txt"

	output:
		R1=GS.remote(GS_PREFIX + "/trimmed/{id}_1.t.fastq.gz"),
		R2=GS.remote(GS_PREFIX + "/trimmed/{id}_2.t.fastq.gz")
	params:
		id="{id}"
	conda: "envs/cutadapt.yaml"
	threads: 4
	shell: "curl https://raw.githubusercontent.com/WatsonLab/GoogleMAGs/master/scripts/ftp_n_trimm.sh | bash -s {params.id} {output.R1} {output.R2}"


rule megahit:
	input:
		R1=GS.remote(GS_PREFIX + "/trimmed/{id}_1.t.fastq.gz"),
		R2=GS.remote(GS_PREFIX + "/trimmed/{id}_2.t.fastq.gz")
	params:
		di=GS.remote(GS_PREFIX + "/megahit/{id}")
	output: 
#		di=GS.remote(GS_PREFIX + "/megahit/{id}/"),
		fa=GS.remote(GS_PREFIX + "/megahit/{id}/final.contigs.fa")
	conda: "envs/megahit.yaml"
	threads: 8
	shell: "mkdir -p {params.di} && megahit --continue --k-list 27,47,67,87 --kmin-1pass -m 0.95 --min-contig-len 1000 -t {threads} -1 {input.R1} -2 {input.R2} -o {params.di}"


rule bwa_index:
	input:  GS.remote(GS_PREFIX + "/megahit/{id}/final.contigs.fa")
	output: 
		ann=GS.remote(GS_PREFIX + "/bwa_indices/{id}.fa.ann"),
		pac=GS.remote(GS_PREFIX + "/bwa_indices/{id}.fa.pac"),
		amb=GS.remote(GS_PREFIX + "/bwa_indices/{id}.fa.amb"),
		bwt=GS.remote(GS_PREFIX + "/bwa_indices/{id}.fa.bwt"),
		sa =GS.remote(GS_PREFIX + "/bwa_indices/{id}.fa.sa")
	params:
		idx=GS.remote(GS_PREFIX + "/bwa_indices/{id}.fa")
	conda: "envs/bwa.yaml"
	threads: 8
	shell:
		'''
		bwa index -p {params.idx} {input}
		'''

rule bwa_mem:
	input:
		R1=GS.remote(GS_PREFIX + "/trimmed/{id}_1.t.fastq.gz"),
		R2=GS.remote(GS_PREFIX + "/trimmed/{id}_2.t.fastq.gz"),
		ann=GS.remote(GS_PREFIX + "/bwa_indices/{id2}.fa.ann"),
		pac=GS.remote(GS_PREFIX + "/bwa_indices/{id2}.fa.pac"),
		amb=GS.remote(GS_PREFIX + "/bwa_indices/{id2}.fa.amb"),
		bwt=GS.remote(GS_PREFIX + "/bwa_indices/{id2}.fa.bwt"),
		sa =GS.remote(GS_PREFIX + "/bwa_indices/{id2}.fa.sa")
	output: 
		bam=GS.remote(GS_PREFIX + "/bam/{id}.{id2}.bam"),
		bai=GS.remote(GS_PREFIX + "/bam/{id}.{id2}.bam.bai"),
		fla=GS.remote(GS_PREFIX + "/bam/{id}.{id2}.bam.flagstat")
	params:
		idx=GS.remote(GS_PREFIX + "/bwa_indices/{id2}.fa")
	conda: "envs/bwa.yaml"
	threads: 8
	shell: 
		'''
		bwa mem -t 8 {params.idx} {input.R1} {input.R2} | samtools sort -@8 -m 500M -o {output.bam} -
		samtools index {output.bam}

		samtools flagstat {output.bam} > {output.fla}
		'''
	
rule coverage:
	input: 
		bam=GS.remote(GS_PREFIX + "/bam/{id}.{id2}.bam"),
		bai=GS.remote(GS_PREFIX + "/bam/{id}.{id2}.bam.bai")
	output:
		cov=GS.remote(GS_PREFIX + "/coverage/{id}.{id2}.txt")
	conda: "envs/metabat2.yaml"
	shell:
		'''
		jgi_summarize_bam_contig_depths --outputDepth {output.cov} {input.bam}
		'''	

