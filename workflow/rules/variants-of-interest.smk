rule mpileup_variants_of_interest:
    """
    Get allele count tables of variants of choice (specified in config file ("IRmutations.tsv"))
    """
    input:
        bam="results/alignments/{sample}.star.bam" if config['pipeline'] == 'parabricks' else "results/alignments/{sample}.hisat2.bam",
        idx="results/alignments/{sample}.star.bam.bai" if config['pipeline'] == 'parabricks' else "results/alignments/{sample}.hisat2.bam.bai",
    output:
        "results/variantAnalysis/variantsOfInterest/counts/{mut_id}/{sample}_allele_counts.tsv",
    conda:
        "../envs/variants.yaml"
    priority: 10
    log:
        "logs/variantsOfInterestMpileup/{sample}_{mut_id}.log",
    params:
        region=lambda wildcards: mutationData[
            mutationData.mutID == wildcards.mut_id
        ].Location.iloc[0],
        ref=config["reference"]["genome"].rstrip(".gz"),
        basedir=workflow.basedir,
    shell:
        """
        samtools mpileup {input.bam} -r {params.region} -f {params.ref} 2> {log} | 
        python {params.basedir}/scripts/base-parser.py > {output} 2>> {log}
        """


rule allele_balance_variants_of_interest:
    """
    R script to take allele count tables from mpileupVOI rule and output .xlsx report for all mutations of interest
    """
    input:
        counts=expand(
            "results/variantAnalysis/variantsOfInterest/counts/{mut_id}/{sample}_allele_counts.tsv",
            sample=samples,
            mut_id=mutationData.mutID,
        ),
        metadata=config["metadata"],
        mutations=config["VariantsOfInterest"]["path"],
    output:
        expand(
            "results/variantAnalysis/variantsOfInterest/csvs/{mut_id}_alleleBalance.csv",
            mut_id=mutationData.mutID,
        ),
        alleleBalance="results/variantAnalysis/variantsOfInterest/alleleBalance.xlsx",
        mean_alleleBalance="results/variantAnalysis/variantsOfInterest/mean_alleleBalance.xlsx",
    conda:
        "../envs/diffexp.yaml"
    priority: 10
    log:
        "logs/variantsOfInterestAlleleBalance.log",
    script:
        "../scripts/variants-of-interest.R"



rule variants_of_interest_notebook:
    """
    Notebook to plot frequencies of Variants of interest
    """
    input:
        nb = f"{workflow.basedir}/notebooks/variants-of-interest.ipynb",
        kernel = "results/.kernel.set",
        muts = expand(
            "results/variantAnalysis/variantsOfInterest/csvs/{mut_id}_alleleBalance.csv",
            mut_id=mutationData.mutID,
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
