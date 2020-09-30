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
		touch("resources/reference/ht.index")
	log:
		"logs/hisat2/build_index.log"
	params:
		ss="--ss {input.splice_sites}",
		exon="--exon {input.exons}",
		prefix="resources/reference/ht2index/"
	threads:8
	wrapper:
		"0.65.0/bio/hisat2/index"
	
rule HISAT2align:
	input:
		reads=["resources/reads/{sample}_1.fastq.gz", "resources/reads/{sample}_2.fastq.gz"],
		idx="resources/reference/ht2index/"
	output:
		pipe("resources/alignments/{sample}.temp.bam")
	log:
		"logs/hisat2/{sample}_align.log"
	params:
		ref=config['ref']['genome'],
		dta="--dta",
		q="-q",
		rgid="--rg-id {sample}",
		rg="--rg SM:{sample}"
	threads:12
	wrapper:
		"0.65.0/bio/hisat2/align"

rule SortBams:
    input:
        "resources/alignments/{sample}.temp.bam"
    output:
        "resources/alignments/{sample}.bam"
    wildcard_constraints:
        sample="^[A-Z]\d$"
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
		metadata="resources/samples.tsv"
	output:
		bamlist="resources/bam.list",
		pops="resources/populations.tsv"
	shell:
		"""
		ls resources/alignments/*bam > {output.bamlist}
		cut -f 4,7 {input.metadata} | tail -n +2 > {output.pops}
		"""

rule VariantCallingFreebayes:
	input:
		bams=expand("resources/alignments/{sample}.bam", sample=samples),
		index=expand("resources/alignments/{sample}.bam.bai", sample=samples),
		ref=config['ref']['genome'],
		samples="resources/bam.list"
	output:
		"results/variants/variants.{chrom}.vcf"
	log:
		"logs/freebayes/{chrom}.log"
	params:
		ploidy="--ploidy " + str(config['ploidy'])
		chrom="-r {chrom}",
		pooled="--pooled-discrete",
		bestn="--use-best-n-alleles 5",
		pops="--populations resources/populations.tsv"
	threads:1
	wrapper:
		"0.65.0/bio/freebayes"

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
        "results/variants/annot.missense.{chrom}.vcf.gz"
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
        vcf=
    conda:
        "../envs/variants.yaml"
    log:
        "logs/allelicdepth/makebed.{chrom}.log"
    shell:
        "vcf2bed < {input.vcf} | cut -f 1-3 > {output} 2> {log}"
