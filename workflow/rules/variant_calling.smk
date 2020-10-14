################ 	Variant calling 	##################
# HISAT2 (Kim et al., 2019)                                              		
# samtools mpileup (Li et al., 2009)                                             
# FreeBayes (Kim et al., 2018)                                                    


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
		"logs/hisat2/splice_sites.log"
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
		directory("resources/reference/ht2index/")
	log:
		"logs/hisat2/index.log"
	conda:
	    "../envs/variants.yaml"
	params:
		ss="--ss {input.splice_sites}",
		exon="--exon {input.exons}",
		prefix="resources/reference/ht2index/"
	threads:8
	shell:
		"hisat2-build -p {threads} --ss {input.splice_sites} --exon {input.exons} {input.fasta} {params.prefix}  2> {log}"
	
rule HISAT2align:
	input:
		reads=["resources/reads/{sample}_1.fastq.gz", "resources/reads/{sample}_2.fastq.gz"],
		idx="resources/reference/ht2index/"
	output:
		temp("temp/{sample}.bam")
	log:
		"logs/hisat2/{sample}_align.log"
	params:
		extra="--dta -q --rg-id {sample} --rg SM:{sample}"
	threads:12
	wrapper:
		"0.65.0/bio/hisat2/align"

rule SortBams:
    input:
        "temp/{sample}.bam"
    output:
        "resources/alignments/{sample}.bam"
    log:
        "logs/samtools/sortbams/{sample}.log"
    wrapper:
        "0.65.0/bio/samtools/sort"


rule IndexBams:
    input:
        "resources/alignments/{sample}.bam"
    output:
        "resources/alignments/{sample}.bam.bai"
    log:
        "logs/samtools/indexbams/{sample}.log"
    wrapper:
        "0.65.0/bio/samtools/index"


rule GenerateParamsFreebayes:
	input:
		ref_idx=config['ref']['genome'],
		metadata=config['samples'],
		bams = expand("resources/alignments/{sample}.bam", sample=samples)
	output:
		bamlist="resources/bam.list",
		pops="resources/populations.tsv"
	shell:
		"""
		ls resources/alignments/*bam > {output.bamlist}
		cut -f 4,7 {input.metadata} | tail -n +2 > {output.pops}
		"""

rule makeRegions:
    input:
        index=lambda wildcards: config['ref']['genome'] + ".fai"
    output:
        expand("resources/regions/genome.{chrom}.region.{i}.bed", chrom=config['chroms'], i = [1,2,3,4,5])
    params:
        chroms=config['chroms']
    script:
        "../scripts/makeRegions.R"


rule VariantCallingFreebayes:
	input:
		bams=expand("resources/alignments/{sample}.bam", sample=samples),
		index=expand("resources/alignments/{sample}.bam.bai", sample=samples),
		ref=config['ref']['genome'],
	        samples="resources/bam.list",
                regions="resources/regions/genome.{chrom}.region.{i}.bed"
	output:
		"results/variants/variants.{chrom}.{i}.vcf"
	log:
		"logs/freebayes/{chrom}.{i}.log"
	params:
		ploidy=config['ploidy'],
		pops="resources/populations.tsv"
	conda:
		"../envs/variants.yaml"
	threads:1
	shell:	"freebayes -f {input.ref} -t {input.regions} --ploidy {params.ploidy} --populations {params.pops} --pooled-discrete --use-best-n-alleles 5 -L {input.samples} > {output} 2> {log}"

rule bcftoolsconcat:
        input:
                expand("results/variants/variants.{chrom}.{i}.vcf", chrom=config['chroms'], i=[1,2,3,4,5])
        output:
                "results/variants/variants.{chrom}.vcf"
        log:
                "logs/bcftools/{chrom}.log"
        conda:
                "../envs/variants.yaml"
        threads:4
        shell:  "bcftools concat -o {output} {input} --threads {threads}"


rule snpEffDbDownload:
    output:
        touch("workflow/scripts/snpeff/db.dl")
    log:
        "logs/snpEff/snpeff_download.log"
    params:
        reference=config['ref']['snpeffdb']
    shell:
        "java -jar workflow/scripts/snpEff/snpEff.jar download {params.reference} 2> {log}"

rule snpEff:
    input:
        calls="results/variants/variants.{chrom}.vcf",
        touch="workflow/scripts/snpeff/db.dl"
    output:
        calls="results/variants/annot.variants.{chrom}.vcf.gz",
    log:
        "logs/snpEff/snpeff_{chrom}.log"
    params:
        db=config['ref']['snpeffdb'],
        prefix="results/variants/annot.variants.{chrom}.vcf"
    shell:
        """
        java -jar workflow/scripts/snpEff/snpEff.jar eff {params.db} {input.calls} > {params.prefix}
        bgzip {params.prefix}
        """

rule MissenseAndQualFilter:
    input:
        vcf="results/variants/annot.variants.{chrom}.vcf.gz"
    output:
        "results/variants/annot.missense.{chrom}.vcf"
    log:
        "logs/snpSift/missense_vcf_{chrom}.log"
    params:
        expression="(ANN[*].EFFECT has 'missense_variant') & (QUAL >= 30)"
    shell:
        """
        java -jar workflow/scripts/snpEff/SnpSift.jar filter "{params.expression}" {input.vcf} > {output} 2> {log}
        """

rule MakeBedOfMissense:
    input:
        vcf="results/variants/annot.missense.{chrom}.vcf",
    output:
        "results/variants/missense.pos.{chrom}.bed",
    conda:
        "../envs/variants.yaml"
    log:
        "logs/allelicdepth/makebed.{chrom}.log"
    shell:
        "vcf2bed < {input.vcf} | cut -f 1-3 > {output} 2> {log}"
