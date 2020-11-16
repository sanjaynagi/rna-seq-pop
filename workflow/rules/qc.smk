#QC rules 
################################          FASTQC				##################################	                                                    					

rule FastQC:
	input:
		"resources/reads/{sample}_{n}.fastq.gz"
	output:
		html="resources/reads/qc/{sample}_{n}_fastqc.html",
		zip="resources/reads/qc/{sample}_{n}_fastqc.zip"
	log:
		"logs/fastqc/{sample}_{n}_QC.log"
	params:
		outdir="--outdir resources/reads/qc"
	wrapper:
		"0.65.0/bio/fastqc"

rule Coverage:
    input:
        "resources/alignments/{sample}.bam"
    output:
        "resources/alignments/coverage/{sample}.mosdepth.summary.txt"
    log:
        "logs/mosdepth/{sample}.log"
    conda:
        "../envs/depth.yaml"
    params:
        prefix = "resources/alignments/coverage/{sample}"
    threads:4
    shell: "mosdepth --threads {threads} --fast-mode --no-per-base {params.prefix} {input}"
