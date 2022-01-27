
rule GenomeIndex:
    """
    Index the reference genome with samtools
    """
    input:
        ref=config["ref"]["genome"],
    output:
        idx=config["ref"]["genome"] + ".fai",
    log:
        "logs/GenomeIndex.log",
    wrapper:
        "v0.69.0/bio/samtools/faidx"


rule HISAT2index:
    """
    Make a HISAT2 index of the reference genome
    """
    input:
        fasta=config["ref"]["genome"]
    output:
        "resources/reference/ht2index/idx.1.ht2",
        touch("resources/reference/ht2index/.complete"),
    log:
        "logs/HISAT2/HISAT2index.log",
    conda:
        "../envs/variants.yaml"
    params:
        prefix=lambda w, output: output[0].split(os.extsep)[0],
    threads: 8
    shell:
        "hisat2-build -p {threads} {input.fasta} {params.prefix}  2> {log}"


rule HISAT2align:
    """
    Align reads to the genome with HISAT2, mark duplicates with samblaster and sort with samtools
    """
    input:
        reads=getFASTQs,
        idx="resources/reference/ht2index/.complete",
    output:
        "results/alignments/{sample}.bam",
    log:
        align="logs/HISAT2/{sample}_align.log",
        sort="logs/samtoolsSort/{sample}.log",
    conda:
        "../envs/variants.yaml"
    params:
        readflags=lambda wildcards: getFASTQs(
            wildcards=wildcards, rules="HISAT2align"
        ),
        extra="--dta -q --rg-id {sample} --rg SM:{sample} --new-summary",
        idx="resources/reference/ht2index/idx",
    threads: 12
    shell:
        """
        hisat2 {params.extra} --threads {threads} -x {params.idx} {params.readflags} 2> {log.align} | 
        samblaster 2> {log.sort} | samtools sort -@{threads} - -o {output} 2>> {log.sort}
        """


rule IndexBams:
    """
    Index bams with samtools
    """
    input:
        "results/alignments/{sample}.bam",
    output:
        "results/alignments/{sample}.bam.bai",
    log:
        "logs/IndexBams/{sample}.log",
    wrapper:
        "0.65.0/bio/samtools/index"


"""
Get a sequence for the number of chunks to break the genome into
"""
chunks = np.arange(1, config["VariantAnalysis"]["chunks"])


rule GenerateFreebayesParams:
    """
    Set up some parameters for freebayes
    """
    input:
        ref_idx=config["ref"]["genome"],
        index=config["ref"]["genome"] + ".fai",
        bams=expand("results/alignments/{sample}.bam", sample=samples),
    output:
        bamlist="resources/bam.list",
        pops="resources/populations.tsv",
        regions=expand(
            "resources/regions/genome.{chrom}.region.{i}.bed",
            chrom=config["chroms"],
            i=chunks,
        ),
    log:
        "logs/GenerateFreebayesParams.log",
    params:
        metadata=config["metadata"],
        chroms=config["chroms"],
        chunks=config["VariantAnalysis"]["chunks"],
    conda:
        "../envs/diffexp.yaml"
    script:
        "../scripts/GenerateFreebayesParams.R"


rule VariantCallingFreebayes:
    """
    Run freebayes on chunks of the genome, splitting the samples by population (strain)
    """
        input:
            bams=expand("results/alignments/{sample}.bam", sample=samples),
            index=expand("results/alignments/{sample}.bam.bai", sample=samples),
            ref=config["ref"]["genome"],
            samples="resources/bam.list",
            regions=ancient("resources/regions/genome.{chrom}.region.{i}.bed"),
        output:
            temp("results/variantAnalysis/vcfs/{chrom}/variants.{i}.vcf"),
        log:
            "logs/VariantCallingFreebayes/{chrom}.{i}.log",
        params:
            ploidy=config["VariantAnalysis"]["ploidy"],
            pops="resources/populations.tsv",
        conda:
            "../envs/variants.yaml"
        threads: 1
        shell:
            "freebayes -f {input.ref} -t {input.regions} --ploidy {params.ploidy} --populations {params.pops} --pooled-discrete --use-best-n-alleles 5 -L {input.samples} > {output} 2> {log}"


rule ConcatVCFs:
    """
    Concatenate VCFs together
    """
    input:
        calls=expand("results/variantAnalysis/vcfs/{{chrom}}/variants.{i}.vcf", i=chunks),
    output:
        temp("results/variantAnalysis/vcfs/variants.{chrom}.vcf"),
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
        "results/variantAnalysis/vcfs/variants.{chrom}.vcf"
    output:
        temp("results/variantAnalysis/vcfs/variants.{chrom}.snps.vcf")
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
        calls="results/variantAnalysis/vcfs/variants.{chrom}.snps.vcf",
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
        "../envs/fstpca.yaml"
    log:
        "logs/ExtractBedVCF.log",
    params:
        chroms=config["chroms"],
    script:
        "../scripts/ExtractBedVCF.py"
