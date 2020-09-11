#QC rules 


################################          FASTQC				##################################	                                                    					

rule fastqc:
	input:
		"data/reads/{sample}_{n}.fastq.gz"
	output:
		"data/reads/QC/{sample}_{n}_fastqc.html"
	log:
		"logs/fastqc/{sample}_{n}_QC.log"
	params:
		out="data/reads/QC"
	shell:
		"""
		fastqc {input} --outdir {params.out} 2> {log}
		"""
