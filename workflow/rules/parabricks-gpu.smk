
rule star_index:
    input:
        fasta=config["reference"]["genome"].rstrip(".gz"),
        gff= config['reference']['gff'],  
    output:
        directory("resources/reference/star_index"),
    threads: 8
    log:
        "logs/star_index/index.log",
    shell:
        """
        STAR \
         --runThreadN {threads} \
         --runMode genomeGenerate \
         --genomeFastaFiles {input.fasta} \
         --sjdbGTFfile {input.gff} \
         --sjdbGTFtagExonParentTranscript Parent \
         --genomeDir {output} 2> {log} \
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
        recal="results/alignments/recal_gpu_{sample}.txt"
    log:
        "logs/parabricks/fq2bam-{sample}.log",
    params:
        wkdir=wkdir
    shell:
        """
        docker run --rm --gpus all --volume {wkdir} --volume {wkdir}/results/alignments/star \
            --workdir {wkdir} \
            nvcr.io/nvidia/clara/clara-parabricks:4.3.0-1 \
            pbrun rna_fq2bam \
            --in-fq {input.reads} \
            --genome-lib-dir resources/reference/star_index\
            --output-dir results/alignments/star \
            --ref {input.ref_fasta} \
            --out-bam {output.bam} \
            --out-recal-file results/alignments/recal_gpu_{wildcards.sample}.txt \
            --read-files-command zcat
        """

rule haplotype_caller:
    """
    Run a GPU-accelerated haplotypeCaller algorithm.
    """
    input:
        bam = "results/alignments/{sample}.star.bam",
        ref_fasta=config["reference"]["genome"].rstrip(".gz"),
        recal="results/alignments/recal_gpu_{sample}.txt",
    output:
        vcf="results/variantAnalysis/vcfs/haplotypecaller/{contig}/{sample}.g.vcf"
    log:
        "logs/parabricks/haplotypecaller-{contig}-{sample}.log",
    params:
        wkdir=wkdir,
        ploidy=config['VariantAnalysis']['ploidy']
    shell:
        """
        docker run --rm --gpus all --volume {wkdir} --volume {wkdir}
            -w {wkdir} \
            nvcr.io/nvidia/clara/clara-parabricks:4.3.0-1 \
            pbrun haplotypecaller \
            --ref {input.ref_fasta} \
            --intervals {wildcards.contig} \
            --in-bam {input.bam} \
            --in-recal-file {input.recal} \
            --out-variants {output.vcf} \
            --ploidy {ploidy}
        """


rule combine_calls:
    input:
        ref=config['reference']['genome'].rstrip(".gz"),
        gvcfs=expand(
            "results/variantAnalysis/vcfs/haplotypecaller/{{contig}}/{sample}.g.vcf", sample=samples
        ),
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
    params:
        extra=config['reference']['genome'].rstrip(".gz"),
    log:
        "logs/gatk/genotypegvcfs.{contig}.log",
    wrapper:
        "v3.13.1/bio/gatk/genotypegvcfs"