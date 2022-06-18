chunks = np.arange(1, config['VariantAnalysis']['chunks'])

rule ConcatVCFs:
    """
    Concatenate VCFs together
    """
    input:
        calls=expand("results/variantAnalysis/vcfs/freebayes/{{chrom}}/variants.{i}.vcf", i=chunks),
    output:
        temp("results/variantAnalysis/vcfs/freebayes/variants.{chrom}.vcf"),
    log:
        "logs/ConcatVCFs/{chrom}.log",
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
        getVCF
    output:
        temp("results/variantAnalysis/vcfs/variants.{chrom}.vcf")
    log:
        "logs/bcftoolsView/{chrom}.log",
    params:
        extra="-v snps",
    wrapper:
        "v1.0.0/bio/bcftools/view"


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
        ref=config["ref"]["snpeffdb"],
        dir="resources/snpEffdb"
    shell:
        "snpEff download {params.ref} -dataDir {params.dir} 2> {log}"


rule snpEff:
    """
    Run snpEff on the VCFs 
    """
    input:
        calls="results/variantAnalysis/vcfs/variants.{chrom}.vcf",
        dl="workflow/scripts/snpEff/db.dl",
    output:
        calls=expand("results/variantAnalysis/vcfs/{dataset}.{{chrom}}.vcf.gz", dataset=config['dataset']),
        csvStats="results/variantAnalysis/vcfs/snpEff.summary.{chrom}.csv",
    log:
        "logs/snpEff/snpEff.{chrom}.log",
    conda:
        "../envs/snpeff.yaml"
    params:
        db=config["ref"]["snpeffdb"],
        prefix=lambda w, output: os.path.splitext(output[0])[0],
        dir="resources/snpEffdb"
    shell:
        """
        snpEff eff {params.db} -dataDir {params.dir} -csvStats {output.csvStats} {input.calls} > {params.prefix} 2> {log}
        bgzip {params.prefix}
        """


rule MissenseAndQualFilter:
    """
    Filter VCFs for missense variants and quality (used later for diffsnps analysis)
    """
    input:
        vcf="results/variantAnalysis/vcfs/{dataset}.{chrom}.vcf.gz",
    output:
        "results/variantAnalysis/vcfs/annot.missense.{chrom}.vcf",
    log:
        "logs/snpSift/missense_vcf_{chrom}.log",
    conda:
        "../envs/snpeff.yaml"
    params:
        expression="(ANN[*].EFFECT has 'missense_variant') & (QUAL >= 30)",
    shell:
        """
        SnpSift filter "{params.expression}" {input.vcf} > {output} 2> {log}
        """


rule ExtractBedVCF:
    """
    Extract SNP positions from the filtered VCFs (used later for diffsnps analysis)
    """
    input:
        vcf=expand(
            "results/variantAnalysis/vcfs/annot.missense.{chrom}.vcf", chrom=config["chroms"]
        ),
    output:
        bed=expand(
            "resources/regions/missense.pos.{chrom}.bed", chrom=config["chroms"]
        ),
    conda:
        "../envs/pythonGenomics.yaml"
    log:
        "logs/ExtractBedVCF.log",
    params:
        chroms=config["chroms"],
    script:
        "../scripts/ExtractBedVCF.py"
