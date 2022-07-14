
chunks = np.arange(1, config['VariantAnalysis']['chunks'])
rule GenerateFreebayesParams:
    input:
        ref_idx = config['ref']['genome'],
        index = config['ref']['genome'] + ".fai",
        bams = expand("results/alignments/{sample}.bam", sample=samples)
    output:
        bamlist = "results/alignments/bam.list",
        pops = "results/alignments/populations.tsv",
        regions = expand("resources/regions/genome.{chrom}.region.{i}.bed", chrom=config['chroms'], i = chunks)
    log:
        "logs/GenerateFreebayesParams.log"
    params:
        metadata = config['metadata'],
        chroms = config['chroms'],
        chunks = config['VariantAnalysis']['chunks']
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
            samples="results/alignments/bam.list",
            pops= "results/alignments/populations.tsv",
            regions=ancient("resources/regions/genome.{chrom}.region.{i}.bed"),
        output:
            temp("results/variantAnalysis/vcfs/freebayes/{chrom}/variants.{i}.vcf"),
        log:
            "logs/VariantCallingFreebayes/{chrom}.{i}.log",
        params:
            ploidy=config["VariantAnalysis"]["ploidy"],
        conda:
            "../envs/variants.yaml"
        threads: 1
        shell:
            "freebayes -f {input.ref} -t {input.regions} --ploidy {params.ploidy} --populations {input.pops} --pooled-discrete --use-best-n-alleles 5 -L {input.samples} > {output} 2> {log}"

rule octopus:
    input:
        reference = config['ref']['genome'],
        bam=expand("results/alignments/{sample}.bam", sample=samples),
        bai=expand("results/alignments/{sample}.bam.bai", sample=samples)
    output:
        vcf=temp("results/variantAnalysis/vcfs/octopus/variants.{chrom}.vcf"),
        vcf_index="results/variantAnalysis/vcfs/octopus/variants.{chrom}.vcf.tbi",
    params:
        ploidy= config['VariantAnalysis']['ploidy'],
        err_model = "PCRF.X10",
        forest = "resources/forests/germline.v0.7.4.forest",
        max_coverage = 1000
    log:
        "logs/octopus/{chrom}.log"
    threads: 16
    shell:
        """
        octopus -R {input.reference} -I {input.bam} -T {wildcards.chrom} --organism-ploidy {params.ploidy} --downsample-above {params.max_coverage} \
         --sequence-error-model {params.err_model} --forest {params.forest} -o {output} --threads {threads} 2> {log}
        """
