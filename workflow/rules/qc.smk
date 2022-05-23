################################          QC           ##################################


rule CheckInputs:
    """
    Check to see that fastq files exist, and reference files are appropriate
    """
    input:
        ref=config["ref"]["genome"],
    output:
        touch("results/.input.check"),
    params:
        metadata=config["metadata"],
        chroms=config["chroms"],
        gffpath=config["ref"]["gff"],
        gene_names=config["ref"]["genes2transcripts"],
        contrasts=config["contrasts"],
        fastq=config["fastq"]["auto"],
        table=config["fastq"]["table"],
        sweeps=config['miscellaneous']["sweeps"]["activate"],
    log:
        "logs/CheckInputs.log",
    conda:
        "../envs/pythonGenomics.yaml"
    priority: 50
    script:
        "../scripts/checkInputs.py"


rule FastQC:
    """
    QC on fastq read data 
    """
    input:
        getFASTQs,
    output:
        html="resources/reads/qc/{sample}_{n}_fastqc.html",
        zip="resources/reads/qc/{sample}_{n}_fastqc.zip",
    log:
        "logs/FastQC/{sample}_{n}_QC.log",
    params:
        outdir="--outdir resources/reads/qc",
    wrapper:
        "0.74.0/bio/fastqc"




rule cutAdapt:
    input:
        getFASTQs,
    output:
        fastq1="resources/reads/trimmed/{sample}_1.fastq.gz",
        fastq2="resources/reads/trimmed/{sample}_2.fastq.gz",
        qc="resources/trimmed/{sample}.qc.txt"
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters="-a AGAGCACACGTCTGAACTCCAGTCAC -g AGATCGGAAGAGCACACGT -A AGAGCACACGTCTGAACTCCAGTCAC -G AGATCGGAAGAGCACACGT",
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        extra="--minimum-length 1 -q 20"
    log:
        "logs/cutadapt/{sample}.log"
    threads: 4 # set desired number of threads here
    wrapper:
        "v0.86.0/bio/cutadapt/pe"





rule BamStats:
    """
    QC alignment statistics
    """
    input:
        bam="results/alignments/{sample}.bam",
        idx="results/alignments/{sample}.bam.bai",
    output:
        stats="results/alignments/bamStats/{sample}.flagstat",
    log:
        "logs/BamStats/{sample}.log",
    wrapper:
        "0.70.0/bio/samtools/flagstat"


rule Coverage:
    """
    Calculate coverage with mosdepth
    """
    input:
        bam="results/alignments/{sample}.bam",
        idx="results/alignments/{sample}.bam.bai",
    output:
        "results/alignments/coverage/{sample}.mosdepth.summary.txt",
    log:
        "logs/Coverage/{sample}.log",
    conda:
        "../envs/depth.yaml"
    params:
        prefix=lambda w, output: output[0].split(os.extsep)[0],
        windowsize=300,
    threads: 4
    shell:
        "mosdepth --threads {threads} --fast-mode --by {params.windowsize} --no-per-base {params.prefix} {input.bam}"


rule vcfStats:
    """
    QC stats of VCF files
    """
    input:
        vcf=expand("results/variantAnalysis/vcfs/{dataset}.{{chrom}}.vcf.gz", dataset=config['dataset']),
    output:
        vcfStats="results/variantAnalysis/vcfs/stats/{chrom}.txt",
    conda:
        "../envs/variants.yaml"
    log:
        "logs/vcfStats/{chrom}.log",
    shell:
        """
        bcftools stats {input} > {output} 2> {log}
        """

rule Qualimap:
    """
    QC of bam files
    """
    input:
        bam="results/alignments/{sample}.bam",
        idx="results/alignments/{sample}.bam.bai",
        gff=config['ref']['gff']
    output:
        "results/alignments/qualimap/{sample}.pdf",
    log:
        "logs/qualimap/{sample}.log",
    conda:
        "../envs/qc.yaml"
    shell:
        "qualimap bamqc -bam {input.bam} -c -gff {input.gff} -outfile {output} 2> {log}"



rule multiQC:
    """
    Integrate QC statistics from other tools into a final .html report
    """
    input:
        expand("resources/reads/qc/{sample}_{n}_fastqc.zip", sample=samples, n=[1, 2]),
        expand("results/variantAnalysis/vcfs/stats/{chrom}.txt", chrom=config["chroms"]),
        expand(
            "results/alignments/coverage/{sample}.mosdepth.summary.txt", sample=samples
        ),
        expand("results/alignments/bamStats/{sample}.flagstat", sample=samples),
        expand("results/quant/{sample}", sample=samples),
    output:
        "results/multiQC.html",
    params:
        "results/ resources/ logs/",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiQC.log",
    wrapper:
        "0.74.0/bio/multiqc"
