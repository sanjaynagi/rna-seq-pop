################ 	Variant calling and IR mutation reports 	##################
# HISAT2 (Kim et al., 2019)                                              		 #
# samtools mpileup (Li et al., 2009)                                             #
# FreeBayes (Kim et al., 2018)                                                    #
##################################################################################

rule HISAT2splicesites:
	input:
		ref = lambda wildcards:config['ref']['genome'],
		gtf = lambda wildcards:config['ref']['gtf']
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
		fasta = lambda wildcards:config['ref']['genome'],
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
		ref=lambda wildcards:config['ref']['genome'],
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

rule mpileupIR:
    input:
        bam="resources/alignments/{sample}.bam",
        index="resources/alignments/{sample}.bam.bai"
    output:
        "results/allele_balance/counts/{sample}_{mut}_allele_counts.tsv"
    log:
        "logs/mpileup/{sample}_{mut}.log"
    params:
        region = lambda wildcards: mutation_data[mutation_data.Name == wildcards.mut].Location.tolist(),
        ref = lambda wildcards:config['ref']['genome']
    shell:
        """
		mkdir -pv results/allele_balance/counts
        samtools mpileup {input.bam} -r {params.region} -f {params.ref} | python2 workflow/scripts/baseParser.py > {output}
        """

rule AlleleBalance:
    input:
        expand("results/allele_balance/counts/{sample}_{mut}_allele_counts.tsv", sample=samples, mut=mutations)
    output:
        "results/allele_balance/allele_balance.xlsx"
    log:
        "logs/allele_balance/Rscript.log"
    shell:
        """
		mkdir -pv results/allele_balance/csvs
		Rscript workflow/scripts/allele_balance.R 2> {log} 
        """

rule GenerateParamsFreebayes:
	input:
		ref_idx=lambda wildcards:config['ref']['genome'],
		metadata="resources/samples.tsv"
	output:
		regions="/regions/freebayes.regions",
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
		ref=lambda wildcards:config['ref']['genome'],
		samples="resources/bam.list"
	output:
		"results/variants/variants.{chrom}.vcf"
	log:
		"logs/freebayes/{chrom}.log"
	params:
		ploidy="--ploidy 10",
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
        reference=config['ref']['snpeffdb'],
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


#### Differential SNP testing ########

rule missenseFilter:
    input:
        vcf="results/variants/annot.variants.{chrom}.vcf.gz"
    output:
        "results/variants/annot.missense.{chrom}.vcf"
    log:
        "logs/snpsift/missense_vcf_{chrom}.log"
    params:
        expression="ANN[*].EFFECT has 'missense_variant'"
    shell:
        """
        java -jar workflow/scripts/snpEff/SnpSift.jar filter "{params.expression}" {input.vcf} > {output} 2> {log}
        """

rule makeBedOfMissense:
    input:
        vcf="results/variants/annot.missense.{chrom}.vcf",
    output:
        "results/variants/missense.pos.{chrom}.bed"
    log:
        "logs/allelicdepth/makebed.{chrom}.log"
    shell:
        """
        vcf2bed < {input.vcf} | cut -f 1-3 > {output} 2> {log}
        """

rule alleleTable:
    input:
        bam="resources/alignments/{sample}.bam",
        bed="results/variants/missense.pos.{chrom}.bed",
        ref=lambda wildcards:config['ref']['genome']
    output:
        "results/variants/alleleTable/{sample}.chr{chrom}.allele.table"
    log:
        "logs/mpileup/{sample}.{chrom}.log"
    params:
        ignore_indels='false',
        baseflt=-5,
        min_alt=3,
        min_af=0
    shell:
        """
        samtools mpileup -f {input.ref} -l {input.bed} {input.bam} | 
        workflow/scripts/mpileup2readcounts/mpileup2readcounts 0 {params.baseflt} {params.ignore_indels} {params.min_alt} {params.min_af} > {output} 2> {log}
        """






