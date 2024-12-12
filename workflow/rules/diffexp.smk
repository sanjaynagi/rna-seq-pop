
rule KallistoIndex:
    """
    Create a kallisto index of the reference transcriptome
    """
    input:
        fasta=config["reference"]["transcriptome"],
    output:
        index="resources/reference/kallisto.idx",
    group:
        "diffexp"
    log:
        "logs/kallisto/index.log",
    wrapper:
        "v1.15.0/bio/kallisto/index"


rule KallistoQuant:
    """
    Pseudo-align reads for each sample to the reference transcriptome.
    Bootstrap to allow for isoform differential expression.
    """
    input:
        fastq=lambda wildcards: getFASTQs(wildcards=wildcards, rules="KallistoQuant"),
        index="resources/reference/kallisto.idx",
    output:
        directory("results/counts/{sample}"),
    group:
        "diffexp"
    log:
        "logs/kallisto/quant_{sample}.log",
    params:
        extra="-b 100" if config['fastq']['paired'] is True else "-b 100 --single -l 75 -s 5",
    threads: 24
    wrapper:
        "v1.15.0/bio/kallisto/quant"



rule DifferentialGeneExpression:
    """
    Perform differential expression analysis at the gene-level with DESeq2
    Produce PCAs, heatmaps, volcano plots
    """
    input:
        metadata=config["metadata"],
        genes2transcripts=config["reference"]["genes2transcripts"],
        counts=expand("results/counts/{sample}", sample=samples),
    output:
        csvs=expand("results/genediff/{comp}.csv", comp=config["contrasts"]),
        xlsx=expand(
            "results/genediff/{dataset}_diffexp.xlsx", dataset=config["dataset"]
        ),
        countStats="results/counts/countStatistics.tsv",
        normCounts="results/counts/normCounts.tsv",
        kallistoSummary = "results/counts/KallistoQuantSummary.tsv",
        corr_heatmap = "results/counts/heatmap_correlations.png",
        pca="results/counts/PCA.pdf",
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
        "../scripts/diffexp-deseq2-genes.R"

rule DifferentialIsoformExpression:
    """
    Perform differential expression analysis at the isoform-level with Sleuth
    Produce volcano plots
    """
    input:
        metadata=config["metadata"],
        genes2transcripts=config["reference"]["genes2transcripts"],
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
        "../scripts/diffexp-sleuth-isoforms.R"


rule GeneSetEnrichment_notebook:
    """
    Perform hypergeometric test GO terms from a gaf file 
    """
    input:
        nb = f"{workflow.basedir}/notebooks/gene-set-enrichment-analysis.ipynb",
        kernel = "results/.kernel.set",
        vcf=expand(
            "results/variantAnalysis/vcfs/{dataset}.{contig}.vcf.gz",
            contig=config["contigs"],
            dataset=config["dataset"],
        ) if config["VariantAnalysis"]["activate"] else [],
        metadata=config["metadata"],
        gaf=config["DifferentialExpression"]["GSEA"]["gaf"],
        DEresults=expand("results/genediff/{comp}.csv", comp=config["contrasts"]),
        Fst="results/variantAnalysis/selection/FstPerGene.tsv" if config["VariantAnalysis"]['selection']["activate"] else [],        
    output:
        nb = "results/notebooks/gene-set-enrichment-analysis.ipynb",
        docs_nb = "docs/rna-seq-pop-results/notebooks/gene-set-enrichment-analysis.ipynb",
        gsea_de = expand(
            "results/gsea/genediff/{comp}.de.tsv",
            comp=config["contrasts"],
        ),
        gsea_fst = expand(
            "results/gsea/fst/{comp}.fst.tsv",
            comp=config["contrasts"],
        ) if config["VariantAnalysis"]['selection']["activate"] else [],
    log:
        "logs/notebooks/gene-set-enrichment-analysis.log",
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        config = configpath,
        selection=config["VariantAnalysis"]['selection']["activate"]
    shell:
        """
        papermill {input.nb} {output.nb} -k pythonGenomics -p go_path {input.gaf} -p metadata_path {input.metadata} -p config_path {params.config} -p selection {params.selection} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """

# rule geneFamilies_notebook:
#     """
#     Summarise gene expression for gene families
#     """
#     input:
#         nb = f"{workflow.basedir}/notebooks/gene-families-heatmap.ipynb",
#         kernel = "results/.kernel.set",
#         genediff=expand("results/genediff/{comp}.csv", comp=config["contrasts"]),
#         normcounts="results/counts/normCounts.tsv",
#         eggnog=config["miscellaneous"]["GeneFamiliesHeatmap"]["eggnog"],
#         pfam=config["miscellaneous"]["GeneFamiliesHeatmap"]["pfam"],
#         metadata=config["metadata"],
#     output:
#         nb = "results/notebooks/gene-families-heatmap.ipynb",
#         docs_nb = "docs/rna-seq-pop-results/notebooks/gene-families-heatmap.ipynb",
#         heatmaps="results/genediff/GeneFamiliesHeatmap.pdf",
#     log:
#         "logs/notebooks/gene-families-heatmap.log",
#     conda:
#         "../envs/pythonGenomics.yaml"
#     params:
#         config_path = configpath
#     shell:
#         """
#         papermill {input.nb} {output.nb} -k pythonGenomics -p metadata_path {input.metadata} -p config_path {params.config_path} -p go_path {input.eggnog} -p normcounts_path {input.normcounts} -p pfam_path {input.pfam} 2> {log}
#         cp {output.nb} {output.docs_nb} 2>> {log}
#         """



rule counts_qc_notebook:
    input:
        nb = f"{workflow.basedir}/notebooks/counts-qc.ipynb",
        kernel = "results/.kernel.set",
        metadata = config["metadata"],
        countStats = "results/counts/countStatistics.tsv",
        kallistoSummary = "results/counts/KallistoQuantSummary.tsv",
        corr_heatmap = "results/counts/heatmap_correlations.png",
    output:
        nb = "results/notebooks/counts-qc.ipynb",
        docs_nb = "docs/rna-seq-pop-results/notebooks/counts-qc.ipynb"
    conda:
        "../envs/pythonGenomics.yaml"
    log:
        "logs/notebooks/counts-qc.log"
    params:
        wd = wkdir
    shell:
        """
        papermill {input.nb} {output.nb} -k pythonGenomics -p wkdir {params.wd} -p metadata_path {input.metadata} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """


rule diffexp_notebook:
    input:
        nb = f"{workflow.basedir}/notebooks/differential-expression.ipynb",
        kernel = "results/.kernel.set",
        xlsx=expand(
            "results/genediff/{dataset}_diffexp.xlsx", dataset=config["dataset"]
        ) if config["DifferentialExpression"]['gene-level']['activate'] else [],
        isoform_xlsx=expand(
            "results/isoformdiff/{dataset}_isoformdiffexp.xlsx",
            dataset=config["dataset"],
        ) if config["DifferentialExpression"]['isoform-level']['activate'] else [],
    output:
        nb = "results/notebooks/differential-expression.ipynb",
        docs_nb = "docs/rna-seq-pop-results/notebooks/differential-expression.ipynb"
    conda:
        "../envs/pythonGenomics.yaml"
    log:
        "logs/notebooks/differential-expression.log"
    params:
        wd = wkdir,
        dataset = config["dataset"],
    shell:
        """
        papermill {input.nb} {output.nb} -k pythonGenomics -p wkdir {params.wd} -p dataset {params.dataset} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """


