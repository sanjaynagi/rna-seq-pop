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
        "../scripts/snp-statistics.py"

rule pca_notebook:
    """
    Perform pca on genotype data and plot 
    """
    input:
        nb = f"{workflow.basedir}/notebooks/principal-components-analysis.ipynb",
        kernel = "results/.kernel.set",
        vcf=expand(
            "results/variantAnalysis/vcfs/{dataset}.{contig}.vcf.gz",
            contig=config["contigs"],
            dataset=config["dataset"],
        ),
        metadata=config["metadata"],
    output:
        nb = "results/notebooks/principal-components-analysis.ipynb",
        docs_nb = "docs/rna-seq-pop-results/notebooks/principal-components-analysis.ipynb",
        PCAfig=expand(
            "results/variantAnalysis/pca/PCA-{contig}-{dataset}.svg",
            contig=config["contigs"],
            dataset=config["dataset"],
        ),
    log:
        "logs/notebooks/pca.log",
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        dataset=config["dataset"],
        config_path = configpath,
        ploidy=config["VariantAnalysis"]["ploidy"],
        missingprop=config["VariantAnalysis"]["pca"]["missingness"],
        qualflt=30,
    shell:
        """
        papermill {input.nb} {output.nb} -k pythonGenomics -p metadata_path {input.metadata} -p dataset {params.dataset} -p config_path {params.config_path} -p ploidy {params.ploidy} -p missingprop {params.missingprop} -p qualflt {params.qualflt} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """

rule geneticDiversity_notebook:
    """
    Calculate population genetic summary statistics on genotype data 
    """
    input:
        nb = f"{workflow.basedir}/notebooks/genetic-diversity.ipynb",
        kernel = "results/.kernel.set",
        vcf=expand(
            "results/variantAnalysis/vcfs/{dataset}.{contig}.vcf.gz",
            contig=config["contigs"],
            dataset=config["dataset"],
        ),
        metadata=config["metadata"],
    output:
        nb = "results/notebooks/genetic-diversity.ipynb",
        docs_nb = "docs/rna-seq-pop-results/notebooks/genetic-diversity.ipynb",
        inbreedingCoef="results/variantAnalysis/diversity/inbreedingCoef.tsv"
        if config["VariantAnalysis"]["ploidy"] > 1
        else [],
        pi="results/variantAnalysis/diversity/SequenceDiversity.tsv",
        theta="results/variantAnalysis/diversity/WattersonsTheta.tsv",
    log:
        "logs/notebooks/genetic-diversity.log",
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        dataset=config["dataset"],
        config_path = configpath,
        ploidy=config["VariantAnalysis"]["ploidy"],
        missingprop=config["VariantAnalysis"]["geneticDiversity"]["missingness"],
        qualflt=30,
    shell:
        """
        papermill {input.nb} {output.nb} -k pythonGenomics -p metadata_path {input.metadata} -p dataset {params.dataset} -p config_path {params.config_path} -p ploidy {params.ploidy} -p missingprop {params.missingprop} -p qualflt {params.qualflt} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """


rule WindowedFstPBS_notebook:
    """
    Calculate population genetic summary statistics on genotype data 
    """
    input:
        nb = f"{workflow.basedir}/notebooks/windowed-selection.ipynb",
        kernel = "results/.kernel.set",
        vcf=expand(
            "results/variantAnalysis/vcfs/{dataset}.{contig}.vcf.gz",
            contig=config["contigs"],
            dataset=config["dataset"],
        ),
        metadata=config["metadata"],
    output:
        nb = "results/notebooks/windowed-selection.ipynb",
        docs_nb = "docs/rna-seq-pop-results/notebooks/windowed-selection.ipynb",
        Fst=expand(
            "results/variantAnalysis/selection/fst/{wsize}/{comp}.Fst.{contig}.svg",
            comp=config["contrasts"],
            contig=config["contigs"],
            wsize=['1000snp_window', '2000snp_window', '5000snp_window'],
        ),
        PBS=(
            expand(
                "results/variantAnalysis/selection/pbs/{wsize}/{pbscomp}.PBS.{contig}.svg",
                pbscomp=config["VariantAnalysis"]["selection"]["population-branch-statistic"]["contrasts"],
                contig=config["contigs"],
                wsize=['1000snp_window', '2000snp_window', '5000snp_window'],
            )
            if config["VariantAnalysis"]["selection"]["population-branch-statistic"]["activate"]
            else []
        ),
    log:
        "logs/notebooks/windowed-selection.log",
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        dataset=config["dataset"],
        config = configpath,
        pbs=config["VariantAnalysis"]["selection"]["population-branch-statistic"]["activate"],
        pbscomps=config["VariantAnalysis"]["selection"]["population-branch-statistic"]["contrasts"],
        ploidy=config["VariantAnalysis"]["ploidy"],
        missingprop=config["VariantAnalysis"]["selection"]["missingness"],
        qualflt=30,
    shell:
        """
        papermill {input.nb} {output.nb} -k pythonGenomics -p config_path {params.config} -p pbs {params.pbs} -p metadata_path {input.metadata} -p dataset {params.dataset} -p ploidy {params.ploidy} -p missingprop {params.missingprop} -p qualflt {params.qualflt} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """


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
        pbs=config["VariantAnalysis"]["selection"]["population-branch-statistic"]["activate"],
        pbscomps=config["VariantAnalysis"]["selection"]["population-branch-statistic"]["contrasts"],
        contigs=config["contigs"],
        ploidy=config["VariantAnalysis"]["ploidy"],
        missingprop=config["VariantAnalysis"]["selection"]["missingness"],
    script:
        "../scripts/per-gene-fst-pbs.py"


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
        "../scripts/ancestry-informative-markers.py"


contig_2l = 'AgamP4_2L' if '2L' not in config['contigs'] else '2L'
contig_2r = 'AgamP4_2R' if '2R' not in config['contigs'] else '2R'

rule Karyotyping:
    input:
        nb = f"{workflow.basedir}/notebooks/karyotype.ipynb",
        kernel = "results/.kernel.set",
        tagsnps = "resources/karyotype_tag_snps.csv",
        vcf2l=expand("results/variantAnalysis/vcfs/{dataset}.{contig}.vcf.gz",
                        contig=contig_2l,
                        dataset=config["dataset"]),
        vcf2r=expand("results/variantAnalysis/vcfs/{dataset}.{contig}.vcf.gz",
                        contig=contig_2r,
                        dataset=config["dataset"]),
    output:
        nb = "results/notebooks/karyotype.ipynb",
        docs_nb = "docs/rna-seq-pop-results/notebooks/karyotype.ipynb",
        tsv = "results/karyotype/karyotypes.tsv",
        png = "results/karyotype/karyotype_heatmap.png",
    log:
        "logs/notebooks/karyotyping.log",
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        metadata=config["metadata"],
        configpath=configpath,
        ploidy=config["VariantAnalysis"]["ploidy"],
        contig_2r=contig_2r,
        contig_2l=contig_2l,
        inversions=config["VariantAnalysis"]["karyotype"]["inversions"],
        dataset=config["dataset"],
    shell:
        """
        papermill {input.nb} {output.nb} -k pythonGenomics -p contig_r {params.contig_2r} -p contig_l {params.contig_2l} -p tag_snp_path {input.tagsnps} -p config_path {params.configpath} -p metadata_path {params.metadata} -p dataset {params.dataset} -p ploidy {params.ploidy} 2> {log}
        cp {output.nb} {output.docs_nb} 2>> {log}
        """
