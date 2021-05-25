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
        "results/allele_balance/counts/{sample}_{mut}_allele_counts.tsv",
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
            "results/allele_balance/counts/{sample}_{mut}_allele_counts.tsv",
            sample=samples,
            mut=mutationData.Name,
        ),
        samples=config["samples"],
        mutations=config["IRmutations"]["path"],
    output:
        expand(
            "results/allele_balance/csvs/{mut}_allele_balance.csv",
            mut=mutationData.Name,
        ),
        allele_balance="results/allele_balance/allele_balance.xlsx",
        mean_allele_balance="results/allele_balance/mean_allele_balance.xlsx",
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
        "results/variants/alleleTables/{sample}.chr{chrom}.allele.table",
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
        samples=config["samples"],
        gff=config["ref"]["gff"],
        geneNames="resources/gene_names.tsv",
        tables=expand(
            "results/variants/alleleTables/{sample}.chr{chrom}.allele.table",
            sample=samples,
            chrom=config["chroms"],
        ),
    output:
        expand(
            "results/variants/diffsnps/{name}.sig.kissDE.tsv", name=config["contrasts"]
        ),
        expand("results/variants/diffsnps/{name}.kissDE.tsv", name=config["contrasts"]),
        expand(
            "results/variants/diffsnps/{name}.normcounts.tsv", name=config["contrasts"]
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


rule StatisticsAndPCA:
    """
    Calculate population genetic summary statistics and PCA on genotype data 
    """
    input:
        vcf=expand(
            "results/variants/vcfs/annot.variants.{chrom}.vcf.gz",
            chrom=config["chroms"],
        ),
        samples=config["samples"],
        gff=config["ref"]["gff"],
    output:
        PCAfig=expand(
            "results/variants/plots/PCA-{chrom}-{dataset}.png",
            chrom=config["chroms"],
            dataset=config["dataset"],
        ),
        SNPdensityFig=expand(
            "results/variants/plots/{dataset}_SNPdensity_{chrom}.png",
            chrom=config["chroms"],
            dataset=config["dataset"],
        ),
        inbreedingCoef="results/variants/stats/inbreedingCoef.tsv",
        inbreedingCoefMean="results/variants/stats/inbreedingCoef.mean.tsv",
        SequenceDiversity="results/variants/stats/SequenceDiversity.tsv",
    log:
        "logs/StatisticsAndPCA.log",
    conda:
        "../envs/fstpca.yaml"
    params:
        DEcontrasts=config["contrasts"],
        dataset=config["dataset"],
        chroms=config["chroms"],
        ploidy=config["VariantCalling"]["ploidy"],
        missingprop=config["pbs"]["missingness"],
        qualflt=30,
    script:
        "../scripts/StatisticsAndPCA.py"


rule WindowedFstPBS:
    """
    Calculate Fst and PBS in windows
    """
    input:
        samples=config["samples"],
        vcf=expand(
            "results/variants/vcfs/annot.variants.{chrom}.vcf.gz",
            chrom=config["chroms"],
        ),
    output:
        Fst=expand(
            "results/variants/plots/fst/{comp}.{chrom}.fst.{wsize}.png",
            comp=config["contrasts"],
            chrom=config["chroms"],
            wsize=config["pbs"]["windownames"],
        ),
        PBS=(
            expand(
                "results/variants/plots/pbs/{pbscomp}.{chrom}.pbs.{wsize}.png",
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


rule PerGeneFstPBS:
    """
    Calculate Fst and PBS for each gene
    """
    input:
        samples=config["samples"],
        gff=config["ref"]["gff"],
        geneNames="resources/gene_names.tsv",
        vcf=expand(
            "results/variants/vcfs/annot.variants.{chrom}.vcf.gz",
            chrom=config["chroms"],
        ),
    output:
        expand("results/variants/{stat}.tsv", stat=windowedStats),
        "results/variants/TajimasD.tsv",
        "results/variants/SequenceDiv.tsv",
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
            "results/variants/vcfs/annot.variants.{chrom}.vcf.gz",
            chrom=config["chroms"],
        ),
        samples=config["samples"],
        aims_zarr_gambcolu=config["AIMs"]["gambcolu"],
        aims_zarr_arab=config["AIMs"]["arab"],
    output:
        AIMs="results/variants/AIMs/AIMs_summary.tsv",
        AIMs_fig="results/variants/AIMs/AIM_fraction_whole_genome.png",
        n_AIMs="results/variants/AIMs/n_AIMS_per_chrom.tsv",
        AIMs_chroms=expand(
            "results/variants/AIMs/AIM_fraction_{chrom}.tsv", chrom=config["chroms"]
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
            lambda wildcards: "results/variants/vcfs/annot.variants.2L.vcf.gz"
            if wildcards.karyo == "2La"
            else "results/variants/vcfs/annot.variants.2R.vcf.gz"
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
        Fst="results/variants/fst.tsv",
        diffsnps=(
            expand(
                "results/variants/diffsnps/{name}.sig.kissDE.tsv",
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
