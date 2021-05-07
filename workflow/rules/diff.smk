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
		xlsx = expand("results/genediff/{dataset}_diffexp.xlsx", dataset=config['dataset']),
		pca = "results/plots/PCA.pdf",
		countstats = "results/quant/count_statistics.tsv",
		normcounts = "results/quant/normcounts.tsv"
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
		xlsx = expand("results/isoformdiff/{dataset}_isoformdiffexp.xlsx", dataset=config['dataset'])
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
		expand("results/genediff/{dataset}_diffexp.xlsx", dataset = config['dataset']),
		expand("results/isoformdiff/{dataset}_isoformdiffexp.xlsx", dataset = config['dataset'])
	output:
		expand("results/genediff/{name}.{direction}.progressive.tsv", name=config['progressiveGenes']['groups'], direction=["up", "down"]),
		expand("results/isoformdiff/{name}.{direction}.progressive.tsv", name=config['progressiveGenes']['groups'], direction=["up", "down"])
	params:
	    samples = config['samples'],
	    comps = config['progressiveGenes']['groups'],
		pval = config['progressiveGenes']['padj_threshold'],
		fc = config['progressiveGenes']['fc_threshold']
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
		gaf = config['GSEA']['gaf'],
		DEresults = expand("results/genediff/{comp}.csv", comp=config['contrasts']),
		diffsnps = expand("results/variants/diffsnps/{comp}.kissDE.tsv", comp=config['contrasts']) if config['diffsnps']['activate'] else [],
		Fst = "results/variants/fst.tsv",
		PBS = "results/variants/pbs.tsv" if config['pbs']['activate'] else [] 
	output:
		expand("results/gsea/genediff/{comp}.DE.{pathway}.tsv", comp = config['contrasts'], pathway=['kegg', 'GO']),
		expand("results/gsea/fst/{comp}.FST.{pathway}.tsv", comp = config['contrasts'], pathway=['kegg', 'GO']),
		expand("results/gsea/diffsnps/{comp}.diffsnps.{pathway}.tsv", comp = config['contrasts'], pathway=['kegg', 'GO']) if config['diffsnps']['activate'] else [],
		expand("results/gsea/pbs/{pbscomp}.PBS.{pathway}.tsv", pbscomp = config['pbs']['contrasts'], pathway=['kegg', 'GO']) if config['pbs']['activate'] else []
	params:
		VariantCalling = config['VariantCalling']['activate'],
		pbs = config['pbs']['activate'],
		pbscomps = config['pbs']['contrasts'],
		diffsnps = config['diffsnps']['activate'],
		replaceStringKegg = config['GSEA']['replaceString'],
		KeggSpeciesID = config['GSEA']['KeggSpeciesID']
	conda:
		"../envs/gsea.yaml"
	log:
		"logs/GeneSetEnrichment.log"
	script:
		"../scripts/GeneSetEnrichment.R"


rule Ag1000gSweepsDE:
	input:
		DEresults = expand("results/genediff/{comp}.csv", comp=config['contrasts']),
		DEcontrasts = "resources/DE.contrast.list"
	output:
		expand("results/genediff/ag1000gSweeps/{comp}_swept.tsv", comp=config['contrasts'])
	log:
		"logs/Ag1000gSweepsDE.log"
	conda:
		"../envs/fstpca.yaml"
	params:
		pval = config['sweeps']['padj_threshold'],
		fc = config['sweeps']['fc_threshold']
	script:
		"../scripts/Ag1000gSweepsDE.py"


rule GeneCategoryContribution:
	input:
		normcounts = "results/quant/normcounts.tsv",
		samples = config['samples'],
	output:
		"results/quant/percentageContributionGeneCategories.tsv"
	log:
		"logs/GeneCategoryContribution.log"
	conda:
		"../envs/diffexp.yaml"
	script:
		"../scripts/GeneCategoryContribution.R"
