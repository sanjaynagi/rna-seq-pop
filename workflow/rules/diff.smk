################################	Differential Gene expression 	################################
# Kallisto (Bray et al., 2016)                                						               #
# DESeq2 (Love et al., 2014)       
# Sleuth (Pimentel et al., 2017)                       	 						               #
####################################################################################################

rule KallistoIndex:
	input:
		fasta = config['ref']['transcriptome']
	output:
		index = "resources/reference/kallisto.idx"
	group:"diffexp"
	log:
		"logs/kallisto/index.log"
	wrapper:
		"0.66.0/bio/kallisto/index"

rule KallistoQuant:
	input:
		fastq = expand("resources/reads/{{sample}}_{n}.fastq.gz", n=[1,2]),
		index = "resources/reference/kallisto.idx"
	output:
		directory("results/quant/{sample}")
	group:"diffexp"
	log:
		"logs/kallisto/quant_{sample}.log"
	params:
		extra = "-b 100"
	threads:24
	wrapper:
		"0.66.0/bio/kallisto/quant"

rule DifferentialGeneExpression:
	input:
		samples = config['samples'],
		gene_names = config['ref']['genenames'],
		DEcontrasts = "resources/DE.contrast.list",
		counts = expand("results/quant/{sample}", sample=samples)
	output:
		csvs = expand("results/genediff/{comp}.csv", comp=config['contrasts']),
		xlsx = "results/genediff/RNA-Seq_diff.xlsx",
		pca = "results/plots/PCA.pdf",
		countstats = "results/quant/count_statistics.tsv"
	group:"diffexp"
	priority: 10
	conda:
		"../envs/diffexp.yaml"
	log:
		"logs/DifferentialGeneExpression.log"
	script:
		"../scripts/DeseqGeneDE.R"

rule DifferentialIsoformExpression:
	input:
		samples = config['samples'],
		gene_names = config['ref']['genenames'],
		DEcontrasts = "resources/DE.contrast.list",
		counts = expand("results/quant/{sample}", sample=samples)
	output:
                csvs = expand("results/isoformdiff/{comp}.csv", comp=config['contrasts']),
		xlsx = "results/isoformdiff/RNA-Seq_isoformdiff.xlsx"
	group:"diffexp"
	priority: 10
	conda:
		"../envs/sleuth.yaml"
	log:
		"logs/DifferentialIsoformExpression.log"
	script:
	    "../scripts/SleuthIsoformsDE.R"


rule progressiveGenesDE:
	input:
		"results/genediff/RNA-Seq_diff.xlsx",
		"results/isoformdiff/RNA-Seq_isoformdiff.xlsx"
	output:
		expand("results/genediff/{name}.{direction}.progressive.tsv", name=config['progressiveGenes']['groups'], direction=["up", "down"]),
		expand("results/isoformdiff/{name}.{direction}.progressive.tsv", name=config['progressiveGenes']['groups'], direction=["up", "down"])
	params:
	    samples = config['samples'],
	    comps = config['progressiveGenes']['groups']
	conda:
	    "../envs/sleuth.yaml"
	log:
	    "logs/progressiveGenesDE.log"
	script:
	    "../scripts/ProgressiveDE.R"


rule GeneSetEnrichment:
	input:
		samples = config['samples'],
		DEcontrasts = "resources/DE.contrast.list",
		gaf = config['ref']['gaf'],
		DEresults = expand("results/genediff/{comp}.csv", comp=config['contrasts']),
		diffsnps = expand("results/variants/diffsnps/{comp}.kissDE.tsv", comp=config['contrasts']),
		Fst = "results/variants/Fst.tsv",
#		PBS = "results/variants/PBS.tsv"
	output:
		expand("results/gsea/genediff/{comp}.DE.csv", comp = config['contrasts']),
		expand("results/gsea/fst/{comp}.FST.csv", comp = config['contrasts']),
		expand("results/gsea/diffsnps/{comp}.diffsnps.csv", comp = config['contrasts']),
#		expand("results/gsea/pbs/{pbscomp}.PBS.csv", pbscomp = config['pbs']['contrasts'])
	params:
		pbs = config['pbs']['activate'],
		pbscomps = config['pbs']['contrasts'],
		replaceStringKegg = config['GSEA']['replaceString'],
		KeggSpeciesID = config['GSEA']['KeggSpeciesID']
	conda:
		"../envs/gsea.yaml"
	log:
		"logs/GSEA/GeneSetEnrichment_.log"
	script:
		"../scripts/GeneSetEnrichment.R"
