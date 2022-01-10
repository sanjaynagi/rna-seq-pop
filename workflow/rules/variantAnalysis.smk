################ Variant Analysis ##################
# scikit-allel (Miles, Harding) 10.5281/zenodo.3935797

rule SNPstatistics:
    """
    Calculate statistics such as no. of SNPs called in exons/introns/genes
    no. of missing SNPs
    """
    input:
        vcf=expand(
            "results/variantAnalysis/vcfs/{dataset}.{chrom}.vcf.gz",
            chrom=config["chroms"], dataset=config['dataset'],
        ),
        metadata=config["metadata"],
        gff=config["ref"]["gff"],
    output:
        snpsPerGenomicFeature = "results/variantAnalysis/SNPstats/snpsPerGenomicFeature.tsv",
        snpsPerGene = "results/variantAnalysis/SNPstats/nSNPsPerGene.tsv",
        SNPdensityFig=expand(
            "results/variantAnalysis/diversity/{dataset}_SNPdensity_{chrom}.png",
            chrom=config["chroms"],
            dataset=config["dataset"],
        ),
    log:
        "logs/SNPstatistics.log"
    conda:
        "../envs/fstpca.yaml"
    params:
        dataset=config['dataset'],
        chroms=config["chroms"],
        ploidy=config["VariantAnalysis"]["ploidy"],
        missingprop=config["VariantAnalysis"]['selection']["missingness"],
        qualflt=30,
    script:
        "../scripts/SNPstatistics.py"


rule PCA:
    """
    Perform pca on genotype data and plot 
    """
    input:
        vcf=expand(
            "results/variantAnalysis/vcfs/{dataset}.{chrom}.vcf.gz",
            chrom=config["chroms"], dataset=config['dataset'],
        ),
        metadata=config["metadata"]
    output:
        PCAfig=expand(
            "results/variantAnalysis/pca/PCA-{chrom}-{dataset}.png",
            chrom=config["chroms"],
            dataset=config["dataset"],
        ),
    log:
        "logs/pca.log"
    conda:
        "../envs/fstpca.yaml"
    params:
        dataset=config["dataset"],
        chroms=config["chroms"],
        ploidy=config["VariantAnalysis"]["ploidy"],
        missingprop=config['VariantAnalysis']['pca']["missingness"],
        qualflt=30,
    script:
        "../scripts/pca.py"


rule SummaryStatistics:
    """
    Calculate population genetic summary statistics and PCA on genotype data 
    """
    input:
        vcf=expand(
            "results/variantAnalysis/vcfs/{dataset}.{chrom}.vcf.gz",
            chrom=config["chroms"], dataset=config['dataset'],
        ),
        metadata=config["metadata"],
    output:
        inbreedingCoef="results/variantAnalysis/diversity/inbreedingCoef.tsv" if config['VariantAnalysis']['ploidy'] > 1 else [],
        pi="results/variantAnalysis/diversity/SequenceDiversity.tsv",
        pi_png="results/variantAnalysis/diversity/piPerChrom.png",
        theta="results/variantAnalysis/diversity/WattersonsTheta.tsv",
        theta_png="results/variantAnalysis/diversity/thetaPerChrom.png"
    log:
        "logs/SummaryStatistics.log",
    conda:
        "../envs/fstpca.yaml"
    params:
        dataset=config["dataset"],
        chroms=config["chroms"],
        ploidy=config["VariantAnalysis"]["ploidy"],
        missingprop=config['VariantAnalysis']['summaryStatistics']["missingness"],
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
            "results/variantAnalysis/vcfs/{dataset}.{chrom}.vcf.gz",
            chrom=config["chroms"], dataset=config['dataset'],
        ),
    output:
        Fst=expand(
            "results/variantAnalysis/selection/fst/Fst.{comp}.{wsize}.{chrom}.png",
            comp=config["contrasts"],
            chrom=config["chroms"],
            wsize=config['VariantAnalysis']['selection']["pbs"]["windownames"],
        ),
        PBS=(
            expand(
                "results/variantAnalysis/selection/pbs/PBS.{pbscomp}.{wsize}.{chrom}.png",
                pbscomp=config['VariantAnalysis']['selection']["pbs"]["contrasts"],
                chrom=config["chroms"],
                wsize=config['VariantAnalysis']['selection']["pbs"]["windownames"],
            )
            if config['VariantAnalysis']['selection']["pbs"]["activate"]
            else []
        ),
    conda:
        "../envs/fstpca.yaml"
    log:
        "logs/WindowedFstPBS.log",
    params:
        dataset=config['dataset'],
        DEcontrasts=config["contrasts"],
        pbs=config['VariantAnalysis']['selection']["pbs"]["activate"],
        pbscomps=config['VariantAnalysis']['selection']["pbs"]["contrasts"],
        chroms=config["chroms"],
        ploidy=config["VariantAnalysis"]["ploidy"],
        missingprop=config['VariantAnalysis']['selection']["missingness"],
        qualflt=30,
        windowsizes=config['VariantAnalysis']['selection']["pbs"]["windowsizes"],
        windowsteps=config['VariantAnalysis']['selection']["pbs"]["windowsteps"],
        windownames=config['VariantAnalysis']['selection']["pbs"]["windownames"],
    script:
        "../scripts/WindowedFstPBS.py"


rule PerGeneFstPBSDxyPi:
    """
    Calculate Fst and PBS for each gene
    """
    input:
        metadata=config["metadata"],
        gff=config["ref"]["gff"],
        geneNames=config['ref']['genes2transcripts'],
        vcf=expand(
            "results/variantAnalysis/vcfs/{dataset}.{chrom}.vcf.gz",
            chrom=config["chroms"], dataset=config['dataset'],
        ),
    output:
        expand("results/variantAnalysis/selection/{stat}PerGene.tsv", stat=windowedStats),
        "results/variantAnalysis/selection/TajimasDPerGene.tsv",
        "results/variantAnalysis/diversity/SequenceDivPerGene.tsv",
        "results/variantAnalysis/diversity/DxyPerGene.tsv"
    conda:
        "../envs/fstpca.yaml"
    log:
        "logs/PerGeneFstPBSDxyPi.log",
    params:
        dataset=config['dataset'],
        DEcontrasts=config["contrasts"],
        pbs=config['VariantAnalysis']['selection']["pbs"]["activate"],
        pbscomps=config['VariantAnalysis']['selection']["pbs"]["contrasts"],
        chroms=config["chroms"],
        ploidy=config["VariantAnalysis"]["ploidy"],
        missingprop=config['VariantAnalysis']['selection']["missingness"],
    script:
        "../scripts/PerGeneFstPBS.py"


rule AncestryInformativeMarkers:
    """
    Calculate the proportion of An.gambiae / An.coluzzii / An.arabiensis ancestry for each sample
    """
    input:
        vcf=expand(
            "results/variantAnalysis/vcfs/{dataset}.{chrom}.vcf.gz",
            chrom=config["chroms"], dataset=config['dataset'],
        ),
        metadata=config["metadata"],
        aims_zarr_gambcolu=config['VariantAnalysis']["AIMs"]["gambcolu"],
        aims_zarr_arab=config['VariantAnalysis']["AIMs"]["arab"],
    output:
        AIMs="results/variantAnalysis/AIMs/AIMs_summary.tsv",
        AIMs_fig="results/variantAnalysis/AIMs/AIM_fraction_whole_genome.png",
        n_AIMs="results/variantAnalysis/AIMs/n_AIMS_per_chrom.tsv",
        AIMs_chroms=expand(
            "results/variantAnalysis/AIMs/AIM_fraction_{chrom}.tsv", chrom=config["chroms"]
        ),
    log:
        "logs/AncestryInformativeMarkers.log",
    conda:
        "../envs/fstpca.yaml"
    params:
        dataset=config['dataset'],
        chroms=config["chroms"],
        ploidy=config["VariantAnalysis"]["ploidy"],
        missingprop=config['VariantAnalysis']["AIMs"]["missingness"],
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
        "../envs/fstpca.yaml"
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
        expand("results/karyotype/{karyo}.{dataset}.karyo.txt", karyo=config['VariantAnalysis']['karyotype']['inversions'], dataset=config['dataset'])
    output:
        "results/karyotype/karyoFreqs.png"
    log:
        "logs/compKaryo/KaryotypePlots.log",
    conda:
        "../envs/fstpca.yaml"
    params:
        metadata=config['metadata'],
        ploidy=config["VariantAnalysis"]["ploidy"],
        inversions=config['VariantAnalysis']['karyotype']['inversions'],
        dataset=config['dataset']
    script:
        "../scripts/KaryoPlots.py"

