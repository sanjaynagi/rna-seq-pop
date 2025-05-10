rule STARindex:
    """
    Make a STAR index of the reference genome.
    """
    input:
        gtfFile=config["reference"]["gff"].rstrip(".gz"),
        fasta=config["reference"]["genome"].rstrip(".gz"),
    output:
        directory("resources/reference/starindex/"),  # Flagging the output as a directory
        touch("resources/reference/starindex/.complete"),  # Completion tracking
    log:
        "logs/STAR/STARindex.log",
    conda:
        "../envs/variants.yaml"
    params:
        prefix="resources/reference/starindex/",
    threads: 8
    shell:
        """
        STAR --runMode genomeGenerate --genomeDir {params.prefix} --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtfFile} --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100 --runThreadN {threads} --outFileNamePrefix {params.prefix} 2> {log}
        """

rule STARalign:
    """
    Align reads to the genome with STAR, generate an unsorted BAM, and then sort with samtools.
    """
    input:
        reads=lambda wildcards: getFASTQs(wildcards=wildcards, rules="STARalign"),
        idx="resources/reference/starindex/"
    output:
        "results/alignments/{sample}.star.bam"
    log:
        align="logs/STAR/{sample}_align.log",              # Log for alignment
        sort="logs/samtoolsSort/{sample}.log"              # Log for sorting
    conda:
        "../envs/variants.yaml"                            # Conda environment
    params:
        readflags=lambda wildcards: " ".join(getFASTQs(wildcards=wildcards, rules="STARalign")),
        extra="--outSAMtype BAM Unsorted",                # STAR output type
        rg=lambda wildcards: f"--outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} PL:ILLUMINA"
    threads: 12
    shell:
        """
        # Generate unsorted BAM with STAR
        STAR --genomeDir {input.idx} --readFilesIn {params.readflags} --readFilesCommand zcat {params.extra} {params.rg} --runThreadN {threads} --outFileNamePrefix results/alignments/{wildcards.sample}. 2> {log.align}

        # Sort the BAM file with samtools, following output file convention
        samtools sort -@{threads} -o {output} results/alignments/{wildcards.sample}.Aligned.out.bam 2>> {log.sort}

        # Remove intermediate unsorted BAM to save space
        rm results/alignments/{wildcards.sample}.Aligned.out.bam
        """

chunks = np.arange(1, config["VariantAnalysis"]["chunks"])

rule GenerateFreebayesParams:
    input:
        ref_idx=config["reference"]["genome"].rstrip(".gz"),
        index=config["reference"]["genome"].rstrip(".gz") + ".fai",
        bams=expand("results/alignments/{sample}.star.bam", sample=samples),
    output:
        bamlist="results/alignments/bam.list",
        pops="results/alignments/populations.tsv",
        regions=expand(
            "results/variantAnalysis/regions/genome.{contig}.region.{i}.bed",
            contig=config["contigs"],
            i=chunks,
        ),
    log:
        "logs/GenerateFreebayesParams.log",
    params:
        aligner='STAR',
        metadata=config["metadata"],
        contigs=config["contigs"],
        chunks=config["VariantAnalysis"]["chunks"],
    conda:
        "../envs/diffexp.yaml"
    script:
        "../scripts/starGenerateFreebayesParams.R"


rule VariantCallingFreebayes:
    """
    Run freebayes on chunks of the genome, splitting the samples by population (strain)
    """
    input:
        bams=expand("results/alignments/{sample}.star.bam", sample=samples),
        index=expand("results/alignments/{sample}.star.bam.bai", sample=samples),
        ref=config["reference"]["genome"].rstrip(".gz"),
        samples="results/alignments/bam.list",
        pops="results/alignments/populations.tsv",
        regions="results/variantAnalysis/regions/genome.{contig}.region.{i}.bed",
    output:
        temp("results/variantAnalysis/vcfs/freebayes/{contig}/variants.{i}.vcf"),
    log:
        "logs/VariantCallingFreebayes/{contig}.{i}.log",
    params:
        ploidy=config["VariantAnalysis"]["ploidy"],
    conda:
        "../envs/variants.yaml"
    threads: 8
    shell:
        """
        freebayes -f {input.ref} -t {input.regions} --ploidy {params.ploidy} --populations {input.pops} --pooled-discrete --use-best-n-alleles 5 --min-alternate-fraction 0.05 -L {input.samples} > {output} 2> {log}
        """


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
