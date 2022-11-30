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


# rule octopus:
#     input:
#         reference=config["reference"]["genome"],
#         bam=expand("results/alignments/{sample}.bam", sample=samples),
#         bai=expand("results/alignments/{sample}.bam.bai", sample=samples),
#     output:
#         vcf=temp("results/variantAnalysis/vcfs/octopus/variants.{contig}.vcf"),
#         vcf_index="results/variantAnalysis/vcfs/octopus/variants.{contig}.vcf.tbi",
#     params:
#         ploidy=config["VariantAnalysis"]["ploidy"],
#         err_model="PCRF.X10",
#         forest="resources/forests/germline.v0.7.4.forest",
#         max_coverage=1000,
#     log:
#         "logs/octopus/{contig}.log",
#     threads: 16
#     conda:
#         "../envs/pythonGenomics.yaml"
#     shell:
#         """
#         octopus -R {input.reference} -I {input.bam} -T {wildcards.contig} --organism-ploidy {params.ploidy} --downsample-above {params.max_coverage} \
#          --sequence-error-model {params.err_model} --forest {params.forest} -o {output} --threads {threads} 2> {log}
#         """
