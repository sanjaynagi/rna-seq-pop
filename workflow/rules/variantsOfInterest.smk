rule mpileupVariantsOfInterest:
    """
    Get allele count tables of variants of choice (specified in config file ("IRmutations.tsv"))
    """
    input:
        bam="results/alignments/{sample}.star.bam" if config['pipeline'] == 'parabricks' else "results/alignments/{sample}.hisat2.bam",
        idx="results/alignments/{sample}.star.bam.bai" if config['pipeline'] == 'parabricks' else "results/alignments/{sample}.hisat2.bam.bai",
    output:
        "results/variantAnalysis/variantsOfInterest/counts/{sample}_{mut}_allele_counts.tsv",
    conda:
        "../envs/variants.yaml"
    priority: 10
    log:
        "logs/variantsOfInterestMpileup/{sample}_{mut}.log",
    params:
        region=lambda wildcards: mutationData[
            mutationData.Name == wildcards.mut
        ].Location.tolist(),
        ref=config["reference"]["genome"].rstrip(".gz"),
        basedir=workflow.basedir,
    shell:
        """
        samtools mpileup {input.bam} -r {params.region} -f {params.ref} 2> {log} | 
        python {params.basedir}/scripts/BaseParser.py > {output} 2>> {log}
        """


rule AlleleBalanceVariantsOfInterest:
    """
    R script to take allele count tables from mpileupVOI rule and output .xlsx report for all mutations of interest
    """
    input:
        counts=expand(
            "results/variantAnalysis/variantsOfInterest/counts/{sample}_{mut}_allele_counts.tsv",
            sample=samples,
            mut=mutationData.Name,
        ),
        metadata=config["metadata"],
        mutations=config["VariantsOfInterest"]["path"],
    output:
        expand(
            "results/variantAnalysis/variantsOfInterest/csvs/{mut}_alleleBalance.csv",
            mut=mutationData.Name,
        ),
        alleleBalance="results/variantAnalysis/variantsOfInterest/alleleBalance.xlsx",
        mean_alleleBalance="results/variantAnalysis/variantsOfInterest/mean_alleleBalance.xlsx",
    conda:
        "../envs/diffexp.yaml"
    priority: 10
    log:
        "logs/variantsOfInterestAlleleBalance.log",
    script:
        "../scripts/VariantsOfInterestAlleleBalance.R"



rule VariantsOfInterest_notebook:
    """
    Notebook to plot frequencies of Variants of interest
    """
    input:
        nb = f"{workflow.basedir}/notebooks/variants-of-interest.ipynb",
        kernel = "results/.kernel.set",
        muts = expand(
            "results/variantAnalysis/variantsOfInterest/csvs/{mut}_alleleBalance.csv",
            mut=mutationData.Name,
        ),
        VariantsOfInterest=config["VariantsOfInterest"]["path"],
    output:
        nb = "results/notebooks/variants-of-interest.ipynb",
        docs_nb = "docs/rna-seq-pop-results/notebooks/variants-of-interest.ipynb",
        perSampleHeatmap="results/variantAnalysis/variantsOfInterest/VOI.heatmapPerSample.svg",
        perTreatmentHeatmap="results/variantAnalysis/variantsOfInterest/VOI.heatmapPerTreatment.svg",
    log:
        "logs/notebooks/variants-of-interest.log",
    conda:
        "../envs/pythonGenomics.yaml"
    priority: 10
    shell:
        """
        papermill {input.nb} {output.nb} -k pythonGenomics -p voi_path {input.VariantsOfInterest} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """
