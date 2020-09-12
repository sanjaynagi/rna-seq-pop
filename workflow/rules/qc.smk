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
		outdir="resources/reads/qc"
	wrapper:
		"0.65.0/bio/fastqc"
