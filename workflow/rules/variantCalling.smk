# This .smk rule file performs the following tasks:

# Extracts exons and splice sites from the reference files (HISAT2)
# Builds HISAT2 index (HISAT2)
# Aligns reads to genome ---> .bam file (HISAT2)
# Sorts .bam, indexes (samtools)
# Generates parameters for freebayes (R)
# Runs freebayes in parallel, splitting the genome into chunks (freebayes)
# concatenates resulting vcfs (bcftools)
# estimates impact of SNPs (snpEff)
# filters VCFs (snpSift)
# extracts bed file of VCF snp-sites (Python)


rule HISAT2splicesites:
	input:
		ref = config['ref']['genome'],
		gtf = config['ref']['gtf']
	output:
		splice_sites="resources/reference/splice-sites.gtf",
		exons="resources/reference/exon-sites.gtf"
	conda:
		"../envs/variants.yaml"
	log:
		"logs/HISAT2/HISAT2splicesites.log"
	shell:
		"""
		hisat2_extract_splice_sites.py {input.gtf} > {output.splice_sites} 2> {log}
		hisat2_extract_exons.py {input.gtf} > {output.exons} 2>> {log}
		"""

rule HISAT2index:
    input:
        fasta = config['ref']['genome'],
        splice_sites="resources/reference/splice-sites.gtf",
        exons="resources/reference/exon-sites.gtf"
    output:
        "resources/reference/ht2index/idx.1.ht2",
        touch("resources/reference/ht2index/.complete")
    log:
        "logs/HISAT2/HISAT2index.log"
    conda:
        "../envs/variants.yaml"
    params:
        ss = "--ss {input.splice_sites}",
        exon = "--exon {input.exons}",
        prefix = lambda w, output: output[0].split(os.extsep)[0]
    threads:8
    shell:
        "hisat2-build -p {threads} --ss {input.splice_sites} --exon {input.exons} {input.fasta} {params.prefix}  2> {log}"
	
rule HISAT2align:
	input:
		reads=["resources/reads/{sample}_1.fastq.gz", "resources/reads/{sample}_2.fastq.gz"],
		idxmarker="resources/reference/ht2index/.complete"             
	output:
		temp("temp/{sample}.bam")
	log:
		"logs/HISAT2/{sample}_HISAT2align.log"
	params:
		extra="--dta -q --rg-id {sample} --rg SM:{sample}",
		idx="resources/reference/ht2index/idx"     
	threads:12
	wrapper:
		"v0.69.0/bio/hisat2/align"

rule SortBams:
    input:
        "temp/{sample}.bam"
    output:
        temp("resources/alignments/{sample}.bam")
    log:
        "logs/SortBams/{sample}.log"
    wrapper:
        "0.65.0/bio/samtools/sort"

rule IndexBams:
    input:
        "resources/alignments/{sample}.bam"
    output:
        "resources/alignments/{sample}.bam.bai"
    log:
        "logs/IndexBams/{sample}.log"
    wrapper:
        "0.65.0/bio/samtools/index"

rule GenomeIndex:
    input:
        ref = config['ref']['genome']
    output:
        idx = config['ref']['genome'] + ".fai"
    log: 
        "logs/GenomeIndex.log"
    wrapper: 
        "v0.69.0/bio/samtools/faidx"

rule markDups:
    input:
        bam = "resources/alignments/{sample}.bam"
    output:
        bam = "resources/alignments/{sample}.marked.bam",
        metrics = "resources/alignments/dedup/{sample}.metrics.txt"
    log:
        "logs/markDups/{sample}.log"
    wrapper:
        "0.72.0/bio/picard/markduplicates"

chunks = np.arange(1, config['chunks'])

rule GenerateFreebayesParams:
    input:
        ref_idx = config['ref']['genome'],
        index = config['ref']['genome'] + ".fai",
        bams = expand("resources/alignments/{sample}.marked.bam", sample=samples)
    output:
        bamlist = "resources/bam.list",
        pops = "resources/populations.tsv",
        regions = expand("resources/regions/genome.{chrom}.region.{i}.bed", chrom=config['chroms'], i = chunks)
    log:
        "logs/GenerateFreebayesParams.log"
    params:
        metadata = config['samples'],
        chroms = config['chroms'],
        chunks = config['chunks']
    conda:
        "../envs/diffexp.yaml"
    script:
        "../scripts/GenerateFreebayesParams.R"

rule VariantCallingFreebayes:
	input:
		bams = expand("resources/alignments/{sample}.marked.bam", sample=samples),
		index = expand("resources/alignments/{sample}.bam.bai", sample=samples),
		ref = config['ref']['genome'],
		samples = ancient("resources/bam.list"),
		regions = ancient("resources/regions/genome.{chrom}.region.{i}.bed")
	output:
		temp("results/variants/vcfs/{chrom}/variants.{i}.vcf")
	log:
		"logs/VariantCallingFreebayes/{chrom}.{i}.log"
	params:
		ploidy = config['ploidy'],
		pops = "resources/populations.tsv"
	conda:
		"../envs/variants.yaml"
	threads:1
	shell:	"freebayes -f {input.ref} -t {input.regions} --ploidy {params.ploidy} --populations {params.pops} --pooled-discrete --use-best-n-alleles 5 -L {input.samples} > {output} 2> {log}"

rule ConcatVCFs:
    input:
        calls = expand("results/variants/vcfs/{{chrom}}/variants.{i}.vcf", i=chunks)
    output:
        "results/variants/vcfs/variants.{chrom}.vcf"
    log:
        "logs/ConcatVCFs/{chrom}.log"
    conda:
        "../envs/variants.yaml"
    threads:4
    shell:  
        "bcftools concat {input.calls} | vcfuniq > {output} 2> {log}"

rule snpEffDbDownload:
    output:
        touch("workflow/scripts/snpEff/db.dl")
    log:
        "logs/snpEff/snpEffDbDownload.log"
    conda:
         "../envs/snpeff.yaml"
    params:
        ref = config['ref']['snpeffdb']
    shell:
        "snpEff download {params.ref} 2> {log}"

rule snpEff:
    input:
        calls = "results/variants/vcfs/variants.{chrom}.vcf",
        dl = "workflow/scripts/snpEff/db.dl"
    output:
        calls = "results/variants/vcfs/annot.variants.{chrom}.vcf.gz",
    log:
        "logs/snpEff/snpEff.{chrom}.log"
    conda:
         "../envs/snpeff.yaml"
    params:
        db = config['ref']['snpeffdb'],
        prefix = lambda w, output: os.path.splitext(output[0])[0]
    shell:
        """
        snpEff eff {params.db} {input.calls} > {params.prefix} 2> {log}
        bgzip {params.prefix}
        """

rule MissenseAndQualFilter:
    input:
        vcf = "results/variants/vcfs/annot.variants.{chrom}.vcf.gz"
    output:
        "results/variants/vcfs/annot.missense.{chrom}.vcf"
    log:
        "logs/snpSift/missense_vcf_{chrom}.log"
    conda:
         "../envs/snpeff.yaml"
    params:
        expression = "(ANN[*].EFFECT has 'missense_variant') & (QUAL >= 30)"
    shell:
        """
        SnpSift filter "{params.expression}" {input.vcf} > {output} 2> {log}
        """

rule ExtractBedVCF:
    input:
        vcf = expand("results/variants/vcfs/annot.missense.{chrom}.vcf", chrom = config['chroms'])
    output:
        bed = expand("resources/regions/missense.pos.{chrom}.bed", chrom = config['chroms'])
    conda:
        "../envs/fstpca.yaml"
    log:
        "logs/ExtractBedVCF.log"
    params:
        chroms = config['chroms']
    script:
        "../scripts/ExtractBedVCF.py"
