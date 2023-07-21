chunks = np.arange(1, config["VariantAnalysis"]["chunks"])

rule ConcatVCFs:
    """
    Concatenate VCFs together
    """
    input:
        calls=expand(
            "results/variantAnalysis/vcfs/freebayes/{{contig}}/variants.{i}.vcf",
            i=chunks,
        ),
    output:
        temp("results/variantAnalysis/vcfs/freebayes/variants.{contig}.vcf"),
    log:
        "logs/ConcatVCFs/{contig}.log",
    conda:
        "../envs/variants.yaml"
    threads: 4
    shell:
        "bcftools concat {input.calls} | vcfuniq > {output} 2> {log}"


rule RestrictToSNPs:
    """"
    Filter out indels
    """
    input:
        "results/variantAnalysis/vcfs/freebayes/variants.{contig}.vcf",
    output:
        temp("results/variantAnalysis/vcfs/variants.{contig}.vcf"),
    log:
        "logs/bcftoolsView/{contig}.log",
    params:
        extra="-v snps",
    wrapper:
        "v1.15.0/bio/bcftools/view"


rule snpEffDbDownload:
    """
    Download the snpEff database for your species
    """
    output:
        touch("workflow/scripts/snpEff/db.dl"),
    log:
        "logs/snpEff/snpEffDbDownload.log",
    conda:
        "../envs/snpeff.yaml"
    params:
        ref=config["reference"]["snpeffdb"],
        dir="resources/snpEffdb",
    shell:
        "snpEff download {params.ref} -dataDir {params.dir} 2> {log}"


rule snpEff:
    """
    Run snpEff on the VCFs 
    """
    input:
        calls="results/variantAnalysis/vcfs/variants.{contig}.vcf",
        dl="workflow/scripts/snpEff/db.dl",
    output:
        calls=expand(
            "results/variantAnalysis/vcfs/{dataset}.{{contig}}.vcf.gz",
            dataset=config["dataset"],
        ),
        csvStats="results/variantAnalysis/vcfs/snpEff.summary.{contig}.csv",
    log:
        "logs/snpEff/snpEff.{contig}.log",
    conda:
        "../envs/snpeff.yaml"
    params:
        db=config["reference"]["snpeffdb"],
        prefix=lambda w, output: os.path.splitext(output[0])[0],
        dir="resources/snpEffdb",
    shell:
        """
        snpEff eff {params.db} -dataDir {params.dir} -csvStats {output.csvStats} {input.calls} > {params.prefix} 2> {log}
        bgzip {params.prefix}
        """

