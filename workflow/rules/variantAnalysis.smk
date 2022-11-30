################ Variant Analysis ##################
# scikit-allel (Miles, Harding) 10.5281/zenodo.3935797

rule SNPstatistics:
    """
    Calculate statistics such as no. of SNPs called in exons/introns/genes
    no. of missing SNPs
    """
    input:
        vcf=expand(
            "results/variantAnalysis/vcfs/{dataset}.{contig}.vcf.gz",
            contig=config["contigs"],
            dataset=config["dataset"],
        ),
        metadata=config["metadata"],
        gff=config["reference"]["gff"],
    output:
        snpsPerGenomicFeature="results/variantAnalysis/SNPstats/snpsPerGenomicFeature.tsv",
        snpsPerGene="results/variantAnalysis/SNPstats/nSNPsPerGene.tsv",
        SNPdensityFig=expand(
            "results/variantAnalysis/diversity/{dataset}_SNPdensity_{contig}.svg",
            contig=config["contigs"],
            dataset=config["dataset"],
        ),
    log:
        "logs/SNPstatistics.log",
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        dataset=config["dataset"],
        contigs=config["contigs"],
        ploidy=config["VariantAnalysis"]["ploidy"],
        missingprop=config["VariantAnalysis"]["selection"]["missingness"],
        qualflt=30,
    script:
        "../scripts/SNPstatistics.py"


rule PCA:
    """
    Perform pca on genotype data and plot 
    """
    input:
        vcf=expand(
            "results/variantAnalysis/vcfs/{dataset}.{contig}.vcf.gz",
            contig=config["contigs"],
            dataset=config["dataset"],
        ),
        metadata=config["metadata"],
    output:
        PCAfig=expand(
            "results/variantAnalysis/pca/PCA-{contig}-{dataset}.svg",
            contig=config["contigs"],
            dataset=config["dataset"],
        ),
    log:
        "logs/pca.log",
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        dataset=config["dataset"],
        contigs=config["contigs"],
        ploidy=config["VariantAnalysis"]["ploidy"],
        missingprop=config["VariantAnalysis"]["pca"]["missingness"],
        qualflt=30,
    script:
        "../scripts/pca.py"


rule SummaryStatistics:
    """
    Calculate population genetic summary statistics and PCA on genotype data 
    """
    input:
        vcf=expand(
            "results/variantAnalysis/vcfs/{dataset}.{contig}.vcf.gz",
            contig=config["contigs"],
            dataset=config["dataset"],
        ),
        metadata=config["metadata"],
    output:
        inbreedingCoef="results/variantAnalysis/diversity/inbreedingCoef.tsv"
        if config["VariantAnalysis"]["ploidy"] > 1
        else [],
        pi="results/variantAnalysis/diversity/SequenceDiversity.tsv",
        theta="results/variantAnalysis/diversity/WattersonsTheta.tsv",
    log:
        "logs/SummaryStatistics.log",
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        dataset=config["dataset"],
        contigs=config["contigs"],
        ploidy=config["VariantAnalysis"]["ploidy"],
        missingprop=config["VariantAnalysis"]["summaryStatistics"]["missingness"],
        qualflt=30,
    script:
        "../scripts/SummaryStats.py"


rule WindowedFstPBS:
    """
    Calculate Fst and PBS in windows
    """
    input:
        metadata=config["metadata"],
        vcf=expand(
            "results/variantAnalysis/vcfs/{dataset}.{contig}.vcf.gz",
            contig=config["contigs"],
            dataset=config["dataset"],
        ),
    output:
        Fst=expand(
            "results/variantAnalysis/selection/fst/{wsize}/{comp}.Fst.{contig}.svg",
            comp=config["contrasts"],
            contig=config["contigs"],
            wsize=['1000snp_window', '2000snp_window', '5000snp_window'],
        ),
        PBS=(
            expand(
                "results/variantAnalysis/selection/pbs/{wsize}/{pbscomp}.PBS.{contig}.svg",
                pbscomp=config["VariantAnalysis"]["selection"]["pbs"]["contrasts"],
                contig=config["contigs"],
                wsize=['1000snp_window', '2000snp_window', '5000snp_window'],
            )
            if config["VariantAnalysis"]["selection"]["pbs"]["activate"]
            else []
        ),
    conda:
        "../envs/pythonGenomics.yaml"
    log:
        "logs/WindowedFstPBS.log",
    params:
        dataset=config["dataset"],
        DEcontrasts=config["contrasts"],
        pbs=config["VariantAnalysis"]["selection"]["pbs"]["activate"],
        pbscomps=config["VariantAnalysis"]["selection"]["pbs"]["contrasts"],
        contigs=config["contigs"],
        ploidy=config["VariantAnalysis"]["ploidy"],
        missingprop=config["VariantAnalysis"]["selection"]["missingness"],
        qualflt=30,
    script:
        "../scripts/WindowedFstPBS.py"


rule PerGeneFstPBSDxyPi:
    """
    Calculate Fst and PBS for each gene
    """
    input:
        metadata=config["metadata"],
        gff=config["reference"]["gff"],
        geneNames=config["reference"]["genes2transcripts"],
        vcf=expand(
            "results/variantAnalysis/vcfs/{dataset}.{contig}.vcf.gz",
            contig=config["contigs"],
            dataset=config["dataset"],
        ),
    output:
        expand(
            "results/variantAnalysis/selection/{stat}PerGene.tsv", stat=windowedStats
        ),
        "results/variantAnalysis/selection/TajimasDPerGene.tsv",
        "results/variantAnalysis/diversity/SequenceDivPerGene.tsv",
        "results/variantAnalysis/diversity/DxyPerGene.tsv",
    conda:
        "../envs/pythonGenomics.yaml"
    log:
        "logs/PerGeneFstPBSDxyPi.log",
    params:
        dataset=config["dataset"],
        DEcontrasts=config["contrasts"],
        pbs=config["VariantAnalysis"]["selection"]["pbs"]["activate"],
        pbscomps=config["VariantAnalysis"]["selection"]["pbs"]["contrasts"],
        contigs=config["contigs"],
        ploidy=config["VariantAnalysis"]["ploidy"],
        missingprop=config["VariantAnalysis"]["selection"]["missingness"],
    script:
        "../scripts/PerGeneFstPBS.py"


rule AncestryInformativeMarkers:
    """
    Calculate the proportion of An.gambiae / An.coluzzii / An.arabiensis ancestry for each sample
    """
    input:
        vcf=expand(
            "results/variantAnalysis/vcfs/{dataset}.{contig}.vcf.gz",
            contig=config["contigs"],
            dataset=config["dataset"],
        ),
        metadata=config["metadata"],
        aims_zarr_gambcolu=config["VariantAnalysis"]["ancestry"]["gambcolu"],
        aims_zarr_arab=config["VariantAnalysis"]["ancestry"]["arab"],
    output:
        AIMs="results/variantAnalysis/ancestry/AIMs_summary.tsv",
        AIMs_fig="results/variantAnalysis/ancestry/AIM_fraction_whole_genome.svg",
        n_AIMs="results/variantAnalysis/ancestry/n_AIMS_per_chrom.tsv",
        AIMs_chroms=expand(
            "results/variantAnalysis/ancestry/AIM_fraction_{contig}.tsv",
            contig=config["contigs"],
        ),
    log:
        "logs/AncestryInformativeMarkers.log",
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        dataset=config["dataset"],
        contigs=config["contigs"],
        ploidy=config["VariantAnalysis"]["ploidy"],
        missingprop=config["VariantAnalysis"]["ancestry"]["missingness"],
        qualflt=30,
    script:
        "../scripts/AncestryInformativeMarkers.py"


rule Karyotype:
    """
    Use compkaryo to determine the "average" karyotype in each sample, based on tagging SNPs. 2la and 2rb only reliable marker sets.
    """
    input:
        vcf=(
            lambda wildcards: "results/variantAnalysis/vcfs/{dataset}.2L.vcf.gz"
            if wildcards.karyo == "2La"
            else "results/variantAnalysis/vcfs/{dataset}.2R.vcf.gz"
        ),
    output:
        "results/karyotype/{karyo}.{dataset}.karyo.txt",
    log:
        "logs/compKaryo/{dataset}.karyo.{karyo}.log",
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        ploidy=config["VariantAnalysis"]["ploidy"],
        basedir=workflow.basedir,
    shell:
        """
        paste <(bcftools query -l {input.vcf}) \
        <(python {params.basedir}/scripts/compkaryo/compkaryo/compkaryo.py {input.vcf} {wildcards.karyo} -p {params.ploidy}) | 
        column -s $'\\t' -t | sort -k 1 > {output}
        """


rule KaryotypePlots:
    """
    Plot karyotype results
    """
    input:
        expand(
            "results/karyotype/{karyo}.{dataset}.karyo.txt",
            karyo=config["VariantAnalysis"]["karyotype"]["inversions"],
            dataset=config["dataset"],
        ),
    output:
        "results/karyotype/karyoFreqs.svg",
    log:
        "logs/compKaryo/KaryotypePlots.log",
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        metadata=config["metadata"],
        ploidy=config["VariantAnalysis"]["ploidy"],
        inversions=config["VariantAnalysis"]["karyotype"]["inversions"],
        dataset=config["dataset"],
    script:
        "../scripts/KaryoPlots.py"
