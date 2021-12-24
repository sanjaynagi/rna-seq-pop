################ Variant Analysis ##################
# scikit-allel (Miles, Harding) 10.5281/zenodo.3935797


rule mpileupIR:
    """
    Get allele count tables of variants of choice (specified in config file ("IRmutations.tsv"))
    """
    input:
        bam="results/alignments/{sample}.bam",
        index="results/alignments/{sample}.bam.bai",
    output:
        "results/alleleBalance/counts/{sample}_{mut}_allele_counts.tsv",
    conda:
        "../envs/variants.yaml"
    priority: 10
    log:
        "logs/mpileupIR/{sample}_{mut}.log",
    params:
        region=lambda wildcards: mutationData[
            mutationData.Name == wildcards.mut
        ].Location.tolist(),
        ref=config["ref"]["genome"],
        basedir=workflow.basedir,
    shell:
        """
        samtools mpileup {input.bam} -r {params.region} -f {params.ref} 2> {log} | 
        python2 {params.basedir}/scripts/BaseParser.py > {output} 2>> {log}
        """


rule AlleleBalanceIR:
    """
    R script to take allele count tables from mpileupIR rule and output .xlsx report for all mutations of interest
    """
    input:
        counts=expand(
            "results/alleleBalance/counts/{sample}_{mut}_allele_counts.tsv",
            sample=samples,
            mut=mutationData.Name,
        ),
        metadata=config["samples"],
        mutations=config["IRmutations"]["path"],
    output:
        expand(
            "results/alleleBalance/csvs/{mut}_alleleBalance.csv",
            mut=mutationData.Name,
        ),
        alleleBalance="results/alleleBalance/alleleBalance.xlsx",
        mean_alleleBalance="results/alleleBalance/mean_alleleBalance.xlsx",
    conda:
        "../envs/diffexp.yaml"
    priority: 10
    log:
        "logs/AlleleBalanceIR.log",
    script:
        "../scripts/MutAlleleBalance.R"


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
        {params.basedir}/scripts/mpileup2readcounts/mpileup2readcounts 0 {params.baseflt} {params.ignore_indels} {params.min_alt} {params.min_af} > {output} 2>> {log}
        """


rule DifferentialSNPs:
    """
    Test to see if any alleles are enriched in one condition versus the other
    """
    input:
        metadata=config["samples"],
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

rule SNPstatistics:
    """
    Calculate statistics such as no. of SNPs called in exons/introns/genes
    no. of missing SNPs
    """
    input:
        vcf=expand(
            "results/variantAnalysis/vcfs/annot.variants.{chrom}.vcf.gz",
            chrom=config["chroms"],
        ),
        metadata=config["samples"],
        gff=config["ref"]["gff"],
    output:
        snpsPerGenomicFeature = "results/variantAnalysis/stats/snpsPerGenomicFeature.tsv",
        SNPdensityFig=expand(
            "results/variantAnalysis/diversity/{dataset}_SNPdensity_{chrom}.png",
            chrom=config["chroms"],
            dataset=config["dataset"],
        ),
    log:
        "logs/SNPstatistics.log"
    params:

    script:
        "../scripts/SNPstatistics.py"


rule PCA:
    """
    Perform pca on genotype data and plot 
    """
    input:
        vcf=expand(
            "results/variantAnalysis/vcfs/annot.variants.{chrom}.vcf.gz",
            chrom=config["chroms"],
        ),
        metadata=config["samples"]
    output:
        PCAfig=expand(
            "results/variantAnalysis/pca/PCA-{chrom}-{dataset}.png",
            chrom=config["chroms"],
            dataset=config["dataset"],
        ),
    log:
        "logs/pca.log"
    params:
        dataset=config["dataset"],
        chroms=config["chroms"],
        ploidy=config["VariantCalling"]["ploidy"],
        missingprop=config["pbs"]["missingness"],
        qualflt=30,
    script:
        "../scripts/pca.py"


rule SummaryStatistics:
    """
    Calculate population genetic summary statistics and PCA on genotype data 
    """
    input:
        vcf=expand(
            "results/variantAnalysis/vcfs/annot.variants.{chrom}.vcf.gz",
            chrom=config["chroms"],
        ),
        metadata=config["samples"],
    output:
        inbreedingCoef="results/variantAnalysis/stats/inbreedingCoef.tsv" if config['VariantCalling']['ploidy'] > 1 else [],
        SequenceDiversity="results/variantAnalysis/stats/SequenceDiversity.tsv",
    log:
        "logs/SummaryStatsAndPCA.log",
    conda:
        "../envs/fstpca.yaml"
    params:
        dataset=config["dataset"],
        chroms=config["chroms"],
        ploidy=config["VariantCalling"]["ploidy"],
        missingprop=config["pbs"]["missingness"],
        qualflt=30,
    script:
        "../scripts/SummaryStats.py"


rule WindowedFstPBS:
    """
    Calculate Fst and PBS in windows
    """
    input:
        metadata=config["samples"],
        vcf=expand(
            "results/variantAnalysis/vcfs/annot.variants.{chrom}.vcf.gz",
            chrom=config["chroms"],
        ),
    output:
        Fst=expand(
            "results/variantAnalysis/selection/fst/{comp}.{chrom}.fst.{wsize}.png",
            comp=config["contrasts"],
            chrom=config["chroms"],
            wsize=config["pbs"]["windownames"],
        ),
        PBS=(
            expand(
                "results/variantAnalysis/selection/pbs/{pbscomp}.{chrom}.pbs.{wsize}.png",
                pbscomp=config["pbs"]["contrasts"],
                chrom=config["chroms"],
                wsize=config["pbs"]["windownames"],
            )
            if config["pbs"]["activate"]
            else []
        ),
    conda:
        "../envs/fstpca.yaml"
    log:
        "logs/WindowedFstPCA.log",
    params:
        DEcontrasts=config["contrasts"],
        pbs=config["pbs"]["activate"],
        pbscomps=config["pbs"]["contrasts"],
        chroms=config["chroms"],
        ploidy=config["VariantCalling"]["ploidy"],
        missingprop=config["pbs"]["missingness"],
        qualflt=30,
        windowsizes=config["pbs"]["windowsizes"],
        windowsteps=config["pbs"]["windowsteps"],
        windownames=config["pbs"]["windownames"],
    script:
        "../scripts/WindowedFstPBS.py"


rule PerGeneFstPBSDxyPi:
    """
    Calculate Fst and PBS for each gene
    """
    input:
        metadata=config["samples"],
        gff=config["ref"]["gff"],
        geneNames=config['ref']['genes2transcripts'],
        vcf=expand(
            "results/variantAnalysis/vcfs/annot.variants.{chrom}.vcf.gz",
            chrom=config["chroms"],
        ),
    output:
        expand("results/variantAnalysis/selection/{stat}PerGene.tsv", stat=windowedStats),
        "results/variantAnalysis/selection/TajimasDPerGene.tsv",
        "results/variantAnalysis/diversity/SequenceDivPerGene.tsv",
        "results/variantAnalysis/diversity/DxyPerGene.tsv"
    conda:
        "../envs/fstpca.yaml"
    log:
        "logs/PerGeneFstPBS.log",
    params:
        DEcontrasts=config["contrasts"],
        pbs=config["pbs"]["activate"],
        pbscomps=config["pbs"]["contrasts"],
        chroms=config["chroms"],
        ploidy=config["VariantCalling"]["ploidy"],
        missingprop=0.8,
    script:
        "../scripts/PerGeneFstPBS.py"


rule AncestryInformativeMarkers:
    """
    Calculate the proportion of An.gambiae / An.coluzzii / An.arabiensis ancestry for each sample
    """
    input:
        vcf=expand(
            "results/variantAnalysis/vcfs/annot.variants.{chrom}.vcf.gz",
            chrom=config["chroms"],
        ),
        metadata=config["samples"],
        aims_zarr_gambcolu=config["AIMs"]["gambcolu"],
        aims_zarr_arab=config["AIMs"]["arab"],
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
        chroms=config["chroms"],
        ploidy=config["VariantCalling"]["ploidy"],
        missingprop=config["AIMs"]["missingness"],
        qualflt=30,
    script:
        "../scripts/AncestryInformativeMarkers.py"


rule Karyotype:
    """
    Use compkaryo to determine the "average" karyotype in each sample, based on tagging SNPs. 2la and 2rb only reliable marker sets.
    """
    input:
        vcf=(
            lambda wildcards: "results/variantAnalysis/vcfs/annot.variants.2L.vcf.gz"
            if wildcards.karyo == "2La"
            else "results/variantAnalysis/vcfs/annot.variants.2R.vcf.gz"
        ),
    output:
        "results/karyotype/{karyo}.karyo.txt",
    log:
        "logs/compKaryo/Karyo_{karyo}.log",
    conda:
        "../envs/fstpca.yaml"
    params:
        ploidy=config["VariantCalling"]["ploidy"],
        basedir=workflow.basedir,
    shell:
        """
        paste <(bcftools query -l {input.vcf}) \
        <(python {params.basedir}/scripts/compkaryo/compkaryo/compkaryo.py {input.vcf} {wildcards.karyo} -p {params.ploidy}) | 
        column -s $'\\t' -t | sort -k 1 > {output}
        """


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
            if config["diffsnps"]["activate"]
            else []
        ),
    output:
        "results/RNA-Seq-full.xlsx",
        expand("results/venn/{name}_DE.Fst.venn.png", name=config["contrasts"]),
    conda:
        "../envs/fstpca.yaml"
    log:
        "logs/VennDiagrams.log",
    params:
        DEcontrasts=config["contrasts"],
        diffsnps=config["diffsnps"]["activate"],
        percentile=0.05,
    script:
        "../scripts/VennDiagrams.py"
