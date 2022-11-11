
rule KallistoIndex:
    """
    Create a kallisto index of the reference transcriptome
    """
    input:
        fasta=config["ref"]["transcriptome"],
    output:
        index="resources/reference/kallisto.idx",
    group:
        "diffexp"
    log:
        "logs/kallisto/index.log",
    wrapper:
        "v1.15.0/bio/kallisto/index"


rule KallistoQuant:
    """
    Pseudo-align reads for each sample to the reference transcriptome.
    Bootstrap to allow for isoform differential expression.
    """
    input:
        fastq=lambda wildcards: getFASTQs(wildcards=wildcards, rules="KallistoQuant"),
        index="resources/reference/kallisto.idx",
    output:
        directory("results/quant/{sample}"),
    group:
        "diffexp"
    log:
        "logs/kallisto/quant_{sample}.log",
    params:
        extra="-b 100" if config['fastq']['paired'] is True else "-b 100 --single",
    threads: 24
    wrapper:
        "v1.15.0/bio/kallisto/quant"


rule GenomeIndex:
    """
    Index the reference genome with samtools
    """
    input:
        config["ref"]["genome"],
    output:
        config["ref"]["genome"] + ".fai",
    log:
        "logs/GenomeIndex.log",
    wrapper:
        "v1.15.0/bio/samtools/faidx"


rule HISAT2index:
    """
    Make a HISAT2 index of the reference genome
    """
    input:
        fasta=config["ref"]["genome"],
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
        reads=lambda wildcards: getFASTQs(
            wildcards=wildcards, rules="HISAT2align_input"
        ),
        idx="resources/reference/ht2index/.complete",
    output:
        "results/alignments/{sample}.bam",
    log:
        align="logs/HISAT2/{sample}_align.log",
        sort="logs/samtoolsSort/{sample}.log",
    conda:
        "../envs/variants.yaml"
    params:
        readflags=lambda wildcards: getFASTQs(wildcards=wildcards, rules="HISAT2align"),
        extra="--dta -q --rg-id {sample} --rg SM:{sample} --rg PL:ILLUMINA --new-summary",
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
        "v1.15.0/bio/samtools/index"


### GATK pre-processing ### removed for now for speed
# rule genome_dict:
#     input:
#         ref=config['ref']['genome'],
#     output:
#         dict=touch("resources/reference/.dict.complete")
#     conda:
#         "../envs/variants.yaml"
#     log:
#         "logs/samtools/create_dict.log",
#     params:
#         output = lambda x: os.path.splitext(config['ref']['genome'])[0] + ".dict",
#     shell:
#         "samtools dict {input} > {params.output} 2> {log} "
# rule splitNCigarReads:
#     input:
#         bam="results/alignments/{sample}.bam",
#         ref=config['ref']['genome'],
#         dict="resources/reference/.dict.complete"
#     output:
#         bam=temp("results/alignments/{sample}.split.bam"),
#     log:
#         "logs/gatk/splitNCIGARreads/{sample}.log",
#     params:
#         extra="",  # optional
#         java_opts="",  # optional
#     resources:
#         mem_mb=4096,
#     wrapper:
#         "v1.4.0/bio/gatk/splitncigarreads"
# rule gatkBaseRecalibrator:
#     input:
#         bam="results/alignments/{sample}.split.bam",
#         ref=config['ref']['genome'],
#         dict="resources/reference/.dict.complete",
#         known="resources/ag3_gaardian.biallelic.vcf.gz"
#     output:
#         recal_table="results/alignments/recal/{sample}.grp",
#     log:
#         "logs/gatk/baserecalibrator/{sample}.log",
#     params:
#         extra="",  # optional
#         java_opts="",  # optional
#     resources:
#         mem_mb=4096,
#     wrapper:
#         "v1.4.0/bio/gatk/baserecalibrator"
# rule gatk_applybqsr:
#     input:
#         bam="results/alignments/{sample}.split.bam",
#         bai="results/alignments/{sample}.split.bam.bai",
#         ref=config['ref']['genome'],
#         dict="resources/reference/.dict.complete",
#         recal_table="results/alignments/recal/{sample}.grp",
#     output:
#         bam="results/alignments/{sample}.split.bq.bam",
#     log:
#         "logs/gatk/gatk_applybqsr/{sample}.log",
#     params:
#         extra="",  # optional
#         java_opts="",  # optional
#     resources:
#         mem_mb=4096,
#     wrapper:
#         "v1.4.0/bio/gatk/applybqsr"
# rule IndexBams2:
#     """
#     Index bams with samtools
#     """
#     input:
#         bam="results/alignments/{sample}.split.bam",
#     output:
#         bai="results/alignments/{sample}.split.bam.bai",
#     log:
#         "logs/IndexBams2/{sample}.log",
#     wrapper:
#         "0.65.0/bio/samtools/index"
# rule IndexBams3:
#     """
#     Index bams with samtools
#     """
#     input:
#         bam="results/alignments/{sample}.split.bq.bam",
#     output:
#         bai="results/alignments/{sample}.split.bq.bam.bai",
#     log:
#         "logs/IndexBams2/{sample}.log",
#     wrapper:
#         "0.65.0/bio/samtools/index"
