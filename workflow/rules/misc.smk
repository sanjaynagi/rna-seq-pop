# rule AlleleTables:
#     """
#     Create allele tables for all missense variants for diffsnps analysis 
#     """
#     input:
#         bam="results/alignments/{sample}.bam",
#         bed="resources/regions/missense.pos.{contig}.bed",
#         ref=config["reference"]["genome"],
#     output:
#         "results/variantAnalysis/alleleTables/{sample}.chr{contig}.allele.table",
#     conda:
#         "../envs/variants.yaml"
#     log:
#         "logs/AlleleTables/{sample}.{contig}.log",
#     params:
#         basedir=workflow.basedir,
#         ignore_indels="false",
#         baseflt=-5,
#         min_alt=3,
#         min_af=0,
#     shell:
#         """
#         samtools mpileup -f {input.ref} -l {input.bed} {input.bam} 2> {log} |
#         {params.basedir}/scripts/mpileup2readcounts/mpileup2readcounts.cc 0 {params.baseflt} {params.ignore_indels} {params.min_alt} {params.min_af} > {output} 2>> {log}
#         """


# rule DifferentialSNPs:
#     """
#     Test to see if any alleles are enriched in one condition versus the other
#     """
#     input:
#         metadata=config["metadata"],
#         gff=config["reference"]["gff"],
#         geneNames=config["reference"]["genes2transcripts"],
#         tables=expand(
#             "results/variantAnalysis/alleleTables/{sample}.chr{contig}.allele.table",
#             sample=samples,
#             contig=config["contigs"],
#         ),
#     output:
#         expand(
#             "results/variantAnalysis/diffsnps/{name}.sig.kissDE.tsv",
#             name=config["contrasts"],
#         ),
#         expand(
#             "results/variantAnalysis/diffsnps/{name}.kissDE.tsv",
#             name=config["contrasts"],
#         ),
#         expand(
#             "results/variantAnalysis/diffsnps/{name}.normcounts.tsv",
#             name=config["contrasts"],
#         ),
#     conda:
#         "../envs/diffsnps.yaml"
#     log:
#         "logs/DifferentialSNPs.log",
#     params:
#         DEcontrasts=config["contrasts"],
#         contigs=config["contigs"],
#         mincounts=100,
#         pval_flt=0.001,  # pvalues already adjusted but way want extra filter for sig file
#     script:
#         "../scripts/DifferentialSNPs.R"


comp_list = get_venn_list()


rule VennDiagrams:
    """
    Find intersection of DE analyses between comparisons and plot
    """
    input:
        DE=expand("results/genediff/{comp}.csv", comp=config["contrasts"]),
    output:
        venn=expand(
            "results/genediff/venn/{comparisons}-{{dir_}}-Venn.png",
            comparisons=comp_list,
        ),
    conda:
        "../envs/venn.yaml"
    log:
        "logs/VennDiagrams_{dir_}.log",
    params:
        comparisons=config["contrasts"],
        padj_threshold=config["DifferentialExpression"]["venn"]["padj_threshold"],
        dataset=config["dataset"],
    script:
        "../scripts/VennDiagrams.py"
