################################	Differential Gene expression 	################################
# Kallisto (Bray et al., 2016)                                						               #
# DESeq2 (Love et al., 2014)       
# Sleuth (Pimentel et al., 2017)                       	 						               #
####################################################################################################

rule transcriptome_index:
	input:
		lambda wildcards:config['ref']['transcriptome']
	output:
		"data/reference/kallisto.idx"
	log:
		"logs/kallisto/index.log"
	shell:
		"""
		kallisto index -i {output} {input} 2> {log}
		"""

rule kallisto:
	input:
		fwd="data/reads/{sample}_1.fastq.gz",\
		rev="data/reads/{sample}_2.fastq.gz",
		idx="data/reference/kallisto.idx"
	output:
		directory("analysis/quant/{sample}")
	threads:12
	log:
		"logs/kallisto/quant_{sample}.log"
	shell:
		"""
		kallisto quant -i {input.idx} -t {threads} -o {output} {input.fwd} {input.rev} -b 100 2> {log}
		"""

rule gene_differential_expression:
	input:
		expand("analysis/quant/{sample}", sample=samples)
	output:
		"analysis/diff/RNA-Seq_diff.xlsx"
	log:
		"logs/DESeq2/geneDE.log"
	shell:
		"""
		mkdir -pv analysis/diff
		Rscript analysis/scripts/kallisto_DE.R 2> {log}
		"""

rule transcript_differential_expression:
	input:
		expand("analysis/quant/{sample}", sample=samples)
	output:
		"analysis/isoformdiff/RNA-Seq_isoformdiff.xlsx"
	log:
		"logs/DESeq2/geneDE.log"
	shell:
		"""
		mkdir -pv analysis/isoformdiff
		Rscript analysis/scripts/sleuth_isoforms_DE.R 2> {log}
		"""