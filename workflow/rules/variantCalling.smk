chunks = np.arange(1, config["VariantAnalysis"]["chunks"])

rule GenerateFreebayesParams:
    input:
        ref_idx=config["reference"]["genome"].rstrip(".gz"),
        index=config["reference"]["genome"].rstrip(".gz") + ".fai",
        bams=expand("results/alignments/{sample}.bam", sample=samples),
    output:
        bamlist="results/alignments/bam.list",
        pops="results/alignments/populations.tsv",
        regions=expand(
            "resources/regions/genome.{contig}.region.{i}.bed",
            contig=config["contigs"],
            i=chunks,
        ),
    log:
        "logs/GenerateFreebayesParams.log",
    params:
        metadata=config["metadata"],
        contigs=config["contigs"],
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
        ref=config["reference"]["genome"].rstrip(".gz"),
        samples=ancient("results/alignments/bam.list"),
        pops=ancient("results/alignments/populations.tsv"),
        regions=ancient("resources/regions/genome.{contig}.region.{i}.bed"),
    output:
        temp("results/variantAnalysis/vcfs/freebayes/{contig}/variants.{i}.vcf"),
    log:
        "logs/VariantCallingFreebayes/{contig}.{i}.log",
    params:
        ploidy=config["VariantAnalysis"]["ploidy"],
    conda:
        "../envs/variants.yaml"
    threads: 1
    shell:
        "freebayes -f {input.ref} -t {input.regions} --ploidy {params.ploidy} --populations {input.pops} --pooled-discrete --use-best-n-alleles 5 -L {input.samples} > {output} 2> {log}"
