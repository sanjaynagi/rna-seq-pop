rule AlleleTables:
    """
    Create allele tables for all missense variants for diffsnps analysis 
    """
    input:
        bam="results/alignments/{sample}.bam",
        bed="resources/regions/missense.pos.{chrom}.bed",
        ref=config["ref"]["genome"],
    output:
        "results/variantAnalysis/alleleTables/{sample}.chr{chrom}.allele.table",
    conda:
        "../envs/variants.yaml"
    log:
        "logs/AlleleTables/{sample}.{chrom}.log",
    params:
        basedir=workflow.basedir,
        ignore_indels="false",
        baseflt=-5,
        min_alt=3,
        min_af=0,
    shell:
        """
        samtools mpileup -f {input.ref} -l {input.bed} {input.bam} 2> {log} |
        {params.basedir}/scripts/mpileup2readcounts/mpileup2readcounts.cc 0 {params.baseflt} {params.ignore_indels} {params.min_alt} {params.min_af} > {output} 2>> {log}
        """


rule DifferentialSNPs:
    """
    Test to see if any alleles are enriched in one condition versus the other
    """
    input:
        metadata=config["metadata"],
        gff=config["ref"]["gff"],
        geneNames=config['ref']['genes2transcripts'],
        tables=expand(
            "results/variantAnalysis/alleleTables/{sample}.chr{chrom}.allele.table",
            sample=samples,
            chrom=config["chroms"],
        ),
    output:
        expand(
            "results/variantAnalysis/diffsnps/{name}.sig.kissDE.tsv", name=config["contrasts"]
        ),
        expand("results/variantAnalysis/diffsnps/{name}.kissDE.tsv", name=config["contrasts"]),
        expand(
            "results/variantAnalysis/diffsnps/{name}.normcounts.tsv", name=config["contrasts"]
        ),
    conda:
        "../envs/diffsnps.yaml"
    log:
        "logs/DifferentialSNPs.log",
    params:
        DEcontrasts=config["contrasts"],
        chroms=config["chroms"],
        mincounts=100,
        pval_flt=0.001,  # pvalues already adjusted but way want extra filter for sig file
    script:
        "../scripts/DifferentialSNPs.R"


rule VennDiagrams:
    """
    Find intersection of DE analyses between comparisons and plot
    Not working May 2021, v0.3.0 
    """
    input:
        DE=expand(
            "results/genediff/{dataset}_diffexp.xlsx", dataset=config["dataset"]
        ),
        Fst="results/variantAnalysis/FstPerGene.tsv",
        diffsnps=(
            expand(
                "results/variantAnalysis/diffsnps/{name}.sig.kissDE.tsv",
                name=config["contrasts"],
            )
            if config['miscellaneous']["diffsnps"]["activate"]
            else []
        ),
    output:
        "results/RNA-Seq-full.xlsx",
        expand("results/venn/{name}_DE.Fst.venn.png", name=config["contrasts"]),
    conda:
        "../envs/pythonGenomics.yaml"
    log:
        "logs/VennDiagrams.log",
    params:
        DEcontrasts=config["contrasts"],
        diffsnps=config['miscellaneous']["diffsnps"]["activate"],
        percentile=0.05,
    script:
        "../scripts/VennDiagrams.py"
