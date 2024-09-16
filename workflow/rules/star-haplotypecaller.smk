
rule star_index:
    input:
        fasta=config["reference"]["genome"].rstrip(".gz"),
        gff= config['reference']['gff'],  
    output:
        index = directory("resources/reference/star_index"),
    threads: 8
    log:
        "logs/star_index/index.log",
    params:
        extra="--sjdbGTFtagExonParentTranscript Parent"
    shell:
        """
        mkdir -p {output.index}
        STAR \
         --runThreadN {threads} \
         --runMode genomeGenerate \
         --genomeFastaFiles {input.fasta} \
         --sjdbGTFfile {input.gff} \
         {params.extra} \
         --genomeDir {output.index} \
         --outFileNamePrefix {output.index}/ \
         > {log} 2>&1
        """

rule fq2bam:
    """
    Clara parabricks: This tool is the equivalent of fq2bam for RNA-Seq samples, receiving inputs in FASTQ format, 
    performing alignment with the splice-aware STAR algorithm, optionally marking of duplicate reads, 
    and outputting an aligned BAM file ready for variant and fusion calling. 
    """
    input:
        reads=getFASTQs,
        ref_fasta=config["reference"]["genome"].rstrip(".gz"),
        star_index= "resources/reference/star_index",
    output:
        bam="results/alignments/{sample}.star.bam",
        bai="results/alignments/{sample}.star.bam.bai"
        # recal="results/alignments/recal_gpu_{sample}.txt"
    log:
        "logs/parabricks/fq2bam-{sample}.log",
    threads: 1
    resources:
        all_gpus = 1
    params:
        wkdir=wkdir
    shell:
        """
        docker run --user 1006:1006 --rm --gpus all --volume {wkdir}:/workdir --workdir /workdir nvcr.io/nvidia/clara/clara-parabricks:4.3.0-1 \
            pbrun rna_fq2bam \
            --in-fq /workdir/{input.reads[0]} /workdir/{input.reads[1]} \
            --genome-lib-dir /workdir/resources/reference/star_index\
            --output-dir /workdir \
            --ref /workdir/{input.ref_fasta} \
            --out-bam /workdir/{output.bam} \
            --read-files-command zcat 2> {log}
        """


rule haplotype_caller:
    """
    Run a GPU-accelerated haplotypeCaller algorithm.
    """
    input:
        bam = "results/alignments/{sample}.star.bam",
        ref_fasta=config["reference"]["genome"].rstrip(".gz"),
        # recal="results/alignments/recal_gpu_{sample}.txt",
    output:
        vcf="results/variantAnalysis/vcfs/haplotypecaller/{contig}/{sample}.g.vcf"
    log:
        "logs/parabricks/haplotypecaller-{contig}-{sample}.log",
    threads: 1
    resources:
        all_gpus = 1
    params:
        wkdir=wkdir,
        ploidy=config['VariantAnalysis']['ploidy']
    shell:
        """
        docker run --user 1006:1006 --rm --gpus all --volume {wkdir}:/workdir -w /workdir nvcr.io/nvidia/clara/clara-parabricks:4.3.0-1 \
            pbrun haplotypecaller \
            --rna \
            --ref /workdir/{input.ref_fasta} \
            --interval {wildcards.contig} \
            --in-bam /workdir/{input.bam} \
            --out-variants /workdir/{output.vcf} \
            --ploidy {params.ploidy} 2> {log}
        """


rule create_dict:
    input:
        ref_fasta=config["reference"]["genome"].rstrip(".gz"),
    output:
        ref_dict=config["reference"]["genome"].rstrip(".gz").rstrip(".fa").rstrip(".fasta") + ".dict",
    log:
        "logs/gatk/create_dict.log",
    params:
        extra="",  # optional: extra arguments for picard.
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024,
    wrapper:
        "v4.3.0/bio/picard/createsequencedictionary"


rule bgzip:
    input:
        vcf="results/variantAnalysis/vcfs/haplotypecaller/{contig}/{sample}.g.vcf",
    output:
        tbi="results/variantAnalysis/vcfs/haplotypecaller/{contig}/{sample}.g.vcf.gz",
    params:
        extra="", # optional
    threads: 1
    log:
        "logs/gatk/bgzip/{contig}.{sample}.log",
    wrapper:
        "v4.3.0-25-g7d3aa9d/bio/bgzip"

rule tabix:
    input:
        vcf="results/variantAnalysis/vcfs/haplotypecaller/{contig}/{sample}.g.vcf.gz",
    output:
        tbi="results/variantAnalysis/vcfs/haplotypecaller/{contig}/{sample}.g.vcf.gz.tbi",
    log:
        "logs/gatk/tabix/{contig}.{sample}.log",
    params:
        "-p vcf",
    wrapper:
        "v4.3.0/bio/tabix/index"

rule combine_calls:
    input:
        ref=config['reference']['genome'].rstrip(".gz"),
        gvcfs=expand(
            "results/variantAnalysis/vcfs/haplotypecaller/{{contig}}/{sample}.g.vcf.gz", sample=samples
        ),
        tbi=expand("results/variantAnalysis/vcfs/haplotypecaller/{{contig}}/{sample}.g.vcf.gz.tbi", sample=samples)
    output:
        gvcf="results/variantAnalysis/vcfs/haplotypecaller/{contig}/variants.g.vcf.gz",
    log:
        log="logs/gatk/combinegvcfs-{contig}.log",
    wrapper:
        "v3.13.1/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref=config['reference']['genome'].rstrip(".gz"),
        gvcf="results/variantAnalysis/vcfs/haplotypecaller/{contig}/variants.g.vcf.gz",
    output:
        vcf=temp("results/variantAnalysis/vcfs/haplotypecaller/variants.{contig}.vcf"),
    log:
        "logs/gatk/genotypegvcfs.{contig}.log",
    wrapper:
        "v3.13.1/bio/gatk/genotypegvcfs"