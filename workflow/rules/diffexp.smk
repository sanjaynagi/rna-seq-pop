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
        genes2transcripts=config['ref']['genes2transcripts'],
        counts=expand("results/quant/{sample}", sample=samples),
    output:
        csvs=expand("results/genediff/{comp}.csv", comp=config["contrasts"]),
        xlsx=expand(
            "results/genediff/{dataset}_diffexp.xlsx", dataset=config["dataset"]
        ),
        pca="results/plots/PCA.pdf",
        countStats="results/quant/countStatistics.tsv",
        normCounts="results/quant/normCounts.tsv",
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
        genes2transcripts=config['ref']['genes2transcripts'],
        counts=expand("results/quant/{sample}", sample=samples),
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
            name=config['DifferentialExpression']["progressiveGenes"]["groups"],
            direction=["up", "down"],
        ),
        expand(
            "results/isoformdiff/{name}.{direction}.progressive.tsv",
            name=config['DifferentialExpression']["progressiveGenes"]["groups"],
            direction=["up", "down"],
        ),
    params:
        metadata=config["metadata"],
        comps=config['DifferentialExpression']["progressiveGenes"]["groups"],
        pval=config['DifferentialExpression']["progressiveGenes"]["padj_threshold"],
        fc=config['DifferentialExpression']["progressiveGenes"]["fc_threshold"],
    conda:
        "../envs/sleuth.yaml"
    log:
        "logs/progressiveGenesDE.log",
    script:
        "../scripts/ProgressiveDE.R"


rule GeneSetEnrichment:
    """
    Perform Gene Set Enrichment analysis with fgsea, on GO terms and KEGG pathways
    """
    input:
        metadata=config["metadata"],
        gaf=config['DifferentialExpression']["GSEA"]["gaf"],
        DEresults=expand("results/genediff/{comp}.csv", comp=config["contrasts"]),
        diffsnps=(
            expand(
                "results/variantAnalysis/diffsnps/{comp}.kissDE.tsv", comp=config["contrasts"]
            )
            if config['miscellaneous']["diffsnps"]["activate"]
            else []
        ),
        Fst="results/variantAnalysis/selection/FstPerGene.tsv",
        PBS="results/variantAnalysis/selection/PbsPerGene.tsv" if config['VariantAnalysis']['selection']["pbs"]["activate"] else [],
    output:
        expand(
            "results/gsea/genediff/{comp}.DE.{pathway}.tsv",
            comp=config["contrasts"],
            pathway=["kegg", "GO"],
        ),
        expand(
            "results/gsea/fst/{comp}.FST.{pathway}.tsv",
            comp=config["contrasts"],
            pathway=["kegg", "GO"],
        ),
        expand(
            "results/gsea/diffsnps/{comp}.diffsnps.{pathway}.tsv",
            comp=config["contrasts"],
            pathway=["kegg", "GO"],
        ) if config['miscellaneous']["diffsnps"]["activate"] else [],
        expand(
            "results/gsea/pbs/{pbscomp}.PBS.{pathway}.tsv",
            pbscomp=config["VariantAnalysis"]['selection']["pbs"]["contrasts"],
            pathway=["kegg", "GO"],
        ) if config["VariantAnalysis"]['selection']["pbs"]["activate"] else [],
    params:
        DEcontrasts=config["contrasts"],
        VariantCalling=config["VariantAnalysis"]["activate"],
        pbs=config["VariantAnalysis"]['selection']["pbs"]["activate"],
        pbscomps=config["VariantAnalysis"]['selection']["pbs"]["contrasts"],
        diffsnps=config['miscellaneous']["diffsnps"]["activate"],
        replaceStringKegg=config['DifferentialExpression']["GSEA"]["replaceString"],
        KeggSpeciesID=config['DifferentialExpression']["GSEA"]["KeggSpeciesID"],
    conda:
        "../envs/gsea.yaml"
    log:
        "logs/GeneSetEnrichment.log",
    script:
        "../scripts/GeneSetEnrichment.R"


rule Ag1000gSweepsDE:
    """
    Find differentially expressed genes that also lie underneath selective sweeps in the Ag1000g
    """
    input:
        DEresults=expand("results/genediff/{comp}.csv", comp=config["contrasts"]),
        signals = "resources/signals.csv"
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
        pval=config['miscellaneous']["sweeps"]["padj_threshold"],
        fc=config['miscellaneous']["sweeps"]["fc_threshold"],
    script:
        "../scripts/Ag1000gSweepsDE.py"


# rule GeneCategoryContribution:
#     """
#     Determine the proportion of read counts that come from P450s, COEs, GSTs, etc.
#     """
#     input:
#         normcounts="results/quant/normcounts.tsv",
#         metadata=config["metadata"],
#     output:
#         "results/quant/percentageContributionGeneCategories.tsv",
#     log:
#         "logs/GeneCategoryContribution.log",
#     conda:
#         "../envs/diffexp.yaml"
#     script:
#         "../scripts/GeneCategoryContribution.R"

rule geneFamilies:
    input:
        genediff = expand("results/genediff/{comp}.csv", comp=config["contrasts"]),
        normcounts = "results/quant/normCounts.tsv",
        eggnog = config['miscellaneous']['GeneFamiliesHeatmap']['eggnog'],
        pfam = config['miscellaneous']['GeneFamiliesHeatmap']['pfam']
    output:
        heatmaps="results/genediff/GeneFamiliesHeatmap.pdf"
    log: 
        "logs/geneFamilies.log"
    params:
        DEcontrasts=config["contrasts"],
    script:
        "../scripts/GeneFamiliesHeatmap.py"
