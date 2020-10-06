################################	Differential Gene expression 	################################
# Kallisto (Bray et al., 2016)                                						               #
# DESeq2 (Love et al., 2014)       
# Sleuth (Pimentel et al., 2017)                       	 						               #
####################################################################################################

rule KallistoIndex:
	input:
		fasta=config['ref']['transcriptome']
	output:
		index="resources/reference/kallisto.idx"
	group:"diffexp"
	log:
		"logs/kallisto/index.log"
	wrapper:
		"0.66.0/bio/kallisto/index"

rule KallistoQuant:
	input:
		fastq=expand("resources/reads/{{sample}}_{n}.fastq.gz", n=[1,2]),
		index="resources/reference/kallisto.idx"
	output:
		directory("results/quant/{sample}")
	group:"diffexp"
	log:
		"logs/kallisto/quant_{sample}.log"
	params:
		extra="-b 100"
	threads:24
	wrapper:
		"0.66.0/bio/kallisto/quant"

rule DifferentialGeneExpression:
	input:
		samples = config['samples'],
		gene_names= config['ref']['genenames'],
		DEcontrasts="resources/DE.contrast.list",
		counts=expand("results/quant/{sample}", sample=samples)
	output:
		"results/genediff/RNA-Seq_diff.xlsx",
		"results/PCA.pdf"
	group:"diffexp"
	priority: 10
	conda:
		"../envs/diffexp.yaml"
	log:
		"logs/DESeq2/geneDE.log"
	script:
		"../scripts/kallistoDE.R"

rule DifferentialIsoformExpression:
	input:
		samples = config['samples'],
		gene_names= config['ref']['genenames'],
		DEcontrasts="resources/DE.contrast.list",
		counts=expand("results/quant/{sample}", sample=samples)
	output:
		"results/isoformdiff/RNA-Seq_isoformdiff.xlsx"
	group:"diffexp"
	priority: 10
	conda:
		"../envs/sleuth.yaml"
	log:
		"logs/sleuth/isoformDE.log"
	script:
	    "../scripts/sleuthIsoformsDE.R"
