rule set_kernel:
    input:
        f'{workflow.basedir}/envs/pythonGenomics.yaml'
    output:
        touch("results/.kernel.set")
    conda: f'{workflow.basedir}/envs/pythonGenomics.yaml'
    log:
        "logs/notebooks/set_kernel.log"
    shell: 
        """
        python -m ipykernel install --user --name=pythonGenomics 2> {log}
        """

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


rule IndexBams:
    """
    Index bams with samtools
    """
    input:
        bam="results/alignments/{sample}.star.bam" if config['pipeline'] == 'parabricks' else "results/alignments/{sample}.hisat2.bam",
    output:
        idx="results/alignments/{sample}.star.bam.bai" if config['pipeline'] == 'parabricks' else "results/alignments/{sample}.hisat2.bam.bai",
    log:
        "logs/IndexBams/{sample}.log",
    wrapper:
        "v1.15.0/bio/samtools/index"


rule RestrictToSNPs:
    """"
    Filter out indels
    """
    input:   
        "results/variantAnalysis/vcfs/haplotypecaller/variants.{contig}.vcf" if config['pipeline'] == 'parabricks' else "results/variantAnalysis/vcfs/freebayes/variants.{contig}.vcf",
    output:
        temp("results/variantAnalysis/vcfs/variants.{contig}.vcf"),
    log:
        "logs/bcftoolsView/{contig}.log",
    params:
        extra="-v snps",
    wrapper:
        "v1.15.0/bio/bcftools/view"
