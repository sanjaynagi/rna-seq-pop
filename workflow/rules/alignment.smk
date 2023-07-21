
rule KallistoIndex:
    """
    Create a kallisto index of the reference transcriptome
    """
    input:
        fasta=config["reference"]["transcriptome"],
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
        directory("results/counts/{sample}"),
    group:
        "diffexp"
    log:
        "logs/kallisto/quant_{sample}.log",
    params:
        extra="-b 100" if config['fastq']['paired'] is True else "-b 100 --single -l 75 -s 5",
    threads: 24
    wrapper:
        "v1.15.0/bio/kallisto/quant"

rule GenomeUnzip:
    """
    Index the reference genome with samtools
    """
    input:
        config["reference"]["genome"],
    output:
        config["reference"]["genome"].rstrip(".gz"),
    log:
        "logs/GenomeUnzip.log",
    shell:
        "gzip -d -c {input} > {output} 2> {log}"

rule GenomeIndex:
    """
    Index the reference genome with samtools
    """
    input:
        config["reference"]["genome"].rstrip(".gz"),
    output:
        config["reference"]["genome"].rstrip(".gz") + ".fai",
    log:
        "logs/GenomeIndex.log",
    wrapper:
        "v1.15.0/bio/samtools/faidx"


rule HISAT2index:
    """
    Make a HISAT2 index of the reference genome
    """
    input:
        fasta=config["reference"]["genome"].rstrip(".gz"),
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
        samblaster="" if config['fastq']['paired'] is True else "--ignoreUnmated"
    threads: 12
    shell:
        """
        hisat2 {params.extra} --threads {threads} -x {params.idx} {params.readflags} 2> {log.align} | 
        samblaster {params.samblaster} 2> {log.sort} | samtools sort -@{threads} - -o {output} 2>> {log.sort}
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
