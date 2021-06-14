################################	Differential Gene expression 	################################
# Kallisto (Bray et al., 2016)                                						               #
# DESeq2 (Love et al., 2014)
# Sleuth (Pimentel et al., 2017)                       	 						               #
####################################################################################################


rule KallistoIndex:
    """
    Create a kallisto index of the reference transcriptome
    """
    input:
        fasta=config["ref"]["transcriptome"],
    output:
        index="resources/reference/kallisto.idx",
    group:
        "diffexp"
    log:
        "logs/kallisto/index.log",
    wrapper:
        "0.66.0/bio/kallisto/index"


rule KallistoQuant:
    """
    Pseudo-align reads for each sample to the reference transcriptome.
    Bootstrap to allow for isoform differential expression.
    """
    input:
        fastq=getFASTQs,
        index="resources/reference/kallisto.idx",
    output:
        directory("results/quant/{sample}"),
    group:
        "diffexp"
    log:
        "logs/kallisto/quant_{sample}.log",
    params:
        extra="-b 100",
    threads: 24
    wrapper:
        "0.66.0/bio/kallisto/quant"


rule DifferentialGeneExpression:
    """
    Perform differential expression analysis at the gene-level with DESeq2
    Produce PCAs, heatmaps, volcano plots
    """
    input:
        metadata=config["samples"],
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
    priority: 10
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
        metadata=config["samples"],
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
            name=config["progressiveGenes"]["groups"],
            direction=["up", "down"],
        ),
        expand(
            "results/isoformdiff/{name}.{direction}.progressive.tsv",
            name=config["progressiveGenes"]["groups"],
            direction=["up", "down"],
        ),
    params:
        metadata=config["samples"],
        comps=config["progressiveGenes"]["groups"],
        pval=config["progressiveGenes"]["padj_threshold"],
        fc=config["progressiveGenes"]["fc_threshold"],
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
        metadata=config["samples"],
        gaf=config["GSEA"]["gaf"],
        DEresults=expand("results/genediff/{comp}.csv", comp=config["contrasts"]),
        diffsnps=(
            expand(
                "results/variants/diffsnps/{comp}.kissDE.tsv", comp=config["contrasts"]
            )
            if config["diffsnps"]["activate"]
            else []
        ),
        Fst="results/variants/fst.tsv",
        PBS="results/variants/pbs.tsv" if config["pbs"]["activate"] else [],
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
        ) if config["diffsnps"]["activate"] else [],
        expand(
            "results/gsea/pbs/{pbscomp}.PBS.{pathway}.tsv",
            pbscomp=config["pbs"]["contrasts"],
            pathway=["kegg", "GO"],
        ) if config["pbs"]["activate"] else [],
    params:
        DEcontrasts=config["contrasts"],
        VariantCalling=config["VariantCalling"]["activate"],
        pbs=config["pbs"]["activate"],
        pbscomps=config["pbs"]["contrasts"],
        diffsnps=config["diffsnps"]["activate"],
        replaceStringKegg=config["GSEA"]["replaceString"],
        KeggSpeciesID=config["GSEA"]["KeggSpeciesID"],
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
    output:
        expand(
            "results/genediff/ag1000gSweeps/{comp}_swept.tsv", comp=config["contrasts"]
        ),
    log:
        "logs/Ag1000gSweepsDE.log",
    conda:
        "../envs/fstpca.yaml"
    params:
        DEcontrasts=config["contrasts"],
        pval=config["sweeps"]["padj_threshold"],
        fc=config["sweeps"]["fc_threshold"],
    script:
        "../scripts/Ag1000gSweepsDE.py"


rule GeneCategoryContribution:
    """
    Determine the proportion of read counts that come from P450s, COEs, GSTs, etc.
    """
    input:
        normcounts="results/quant/normcounts.tsv",
        metadata=config["samples"],
    output:
        "results/quant/percentageContributionGeneCategories.tsv",
    log:
        "logs/GeneCategoryContribution.log",
    conda:
        "../envs/diffexp.yaml"
    script:
        "../scripts/GeneCategoryContribution.R"
