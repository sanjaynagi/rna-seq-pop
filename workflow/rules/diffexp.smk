################################	Differential Gene expression 	################################
# Kallisto (Bray et al., 2016)                                						               #
# DESeq2 (Love et al., 2014)
# Sleuth (Pimentel et al., 2017)                       	 						               #
####################################################################################################


rule DifferentialGeneExpression:
    """
    Perform differential expression analysis at the gene-level with DESeq2
    Produce PCAs, heatmaps, volcano plots
    """
    input:
        metadata=config["metadata"],
        genes2transcripts=config["ref"]["genes2transcripts"],
        counts=expand("results/counts/{sample}", sample=samples),
    output:
        csvs=expand("results/genediff/{comp}.csv", comp=config["contrasts"]),
        xlsx=expand(
            "results/genediff/{dataset}_diffexp.xlsx", dataset=config["dataset"]
        ),
        pca="results/counts/PCA.pdf",
        countStats="results/counts/countStatistics.tsv",
        normCounts="results/counts/normCounts.tsv",
    group:
        "diffexp"
    priority: 20
    conda:
        "../envs/diffexp.yaml"
    log:
        "logs/DifferentialGeneExpression.log",
    params:
        DEcontrasts=config["contrasts"],
    script:
        "../scripts/DeseqGeneDE.R"


rule DifferentialIsoformExpression:
    """
    Perform differential expression analysis at the isoform-level with Sleuth
    Produce volcano plots
    """
    input:
        metadata=config["metadata"],
        genes2transcripts=config["ref"]["genes2transcripts"],
        counts=expand("results/counts/{sample}", sample=samples),
    output:
        csvs=expand("results/isoformdiff/{comp}.csv", comp=config["contrasts"]),
        xlsx=expand(
            "results/isoformdiff/{dataset}_isoformdiffexp.xlsx",
            dataset=config["dataset"],
        ),
    group:
        "diffexp"
    priority: 10
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/DifferentialIsoformExpression.log",
    params:
        DEcontrasts=config["contrasts"],
    script:
        "../scripts/SleuthIsoformsDE.R"


rule progressiveGenesDE:
    """
    Determine if any genes are up/down regulated in same direction in multiple comparisons
    """
    input:
        expand("results/genediff/{dataset}_diffexp.xlsx", dataset=config["dataset"]),
        expand(
            "results/isoformdiff/{dataset}_isoformdiffexp.xlsx",
            dataset=config["dataset"],
        ),
    output:
        expand(
            "results/genediff/{name}.{direction}.progressive.tsv",
            name=config["DifferentialExpression"]["progressiveGenes"]["groups"],
            direction=["up", "down"],
        ),
        expand(
            "results/isoformdiff/{name}.{direction}.progressive.tsv",
            name=config["DifferentialExpression"]["progressiveGenes"]["groups"],
            direction=["up", "down"],
        ),
    params:
        metadata=config["metadata"],
        comps=config["DifferentialExpression"]["progressiveGenes"]["groups"],
        pval=config["DifferentialExpression"]["progressiveGenes"]["padj_threshold"],
        fc=config["DifferentialExpression"]["progressiveGenes"]["fc_threshold"],
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/progressiveGenesDE.log",
    script:
        "../scripts/ProgressiveDE.R"


rule GeneSetEnrichment:
    """
    Perform hypergeometric test GO terms from a gaf file 
    """
    input:
        metadata=config["metadata"],
        gaf=config["DifferentialExpression"]["GSEA"]["gaf"],
        DEresults=expand("results/genediff/{comp}.csv", comp=config["contrasts"]),
        Fst="results/variantAnalysis/selection/FstPerGene.tsv" if config["VariantAnalysis"]['selection']["activate"] else [],
    output:
        expand(
            "results/gsea/genediff/{comp}.de.tsv",
            comp=config["contrasts"],
        ),
        expand(
            "results/gsea/fst/{comp}.fst.tsv",
            comp=config["contrasts"],
        ) if config["VariantAnalysis"]['selection']["activate"] else [],
    params:
        DEcontrasts=config["contrasts"],
        selection=config["VariantAnalysis"]['selection']["activate"]
    conda:
        "../envs/pythonGenomics.yaml"
    log:
        "logs/GeneSetEnrichment.log",
    script:
        "../scripts/GeneSetEnrichment.py"


rule Ag1000gSweepsDE:
    """
    Find differentially expressed genes that also lie underneath selective sweeps in the Ag1000g
    """
    input:
        DEresults=expand("results/genediff/{comp}.csv", comp=config["contrasts"]),
        signals="resources/signals.csv",
    output:
        expand(
            "results/genediff/ag1000gSweeps/{comp}_swept.tsv", comp=config["contrasts"]
        ),
    log:
        "logs/Ag1000gSweepsDE.log",
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        DEcontrasts=config["contrasts"],
        pval=config["miscellaneous"]["sweeps"]["padj_threshold"],
        fc=config["miscellaneous"]["sweeps"]["fc_threshold"],
    script:
        "../scripts/Ag1000gSweepsDE.py"


rule geneFamilies:
    input:
        genediff=expand("results/genediff/{comp}.csv", comp=config["contrasts"]),
        normcounts="results/counts/normCounts.tsv",
        eggnog=config["miscellaneous"]["GeneFamiliesHeatmap"]["eggnog"],
        pfam=config["miscellaneous"]["GeneFamiliesHeatmap"]["pfam"],
    output:
        heatmaps="results/genediff/GeneFamiliesHeatmap.pdf",
    conda:
        "../envs/pythonGenomics.yaml"
    log:
        "logs/geneFamilies.log",
    params:
        DEcontrasts=config["contrasts"],
    script:
        "../scripts/GeneFamiliesHeatmap.py"
