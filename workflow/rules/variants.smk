################ 	Variant calling and IR mutation reports 	##################
# HISAT2 (Kim et al., 2019)                                              		 #
# samtools mpileup (Li et al., 2009)                                             #
# FreeBayes (Kim et al., 2018)                                                    #
##################################################################################

rule genome_index:
	input:
		ref = lambda wildcards:config['ref']['genome'],
		gtf = lambda wildcards:config['ref']['gtf']
	output:
		touch("data/reference/ht.index")
	log:
		"logs/hisat2/build_index.log"
	shell:
		"""
		hisat2_extract_splice_sites.py {input.gtf} > data/reference/splice-sites.gtf 2> {log}
		hisat2_extract_exons.py {input.gtf} > data/reference/exon-sites.gtf 2>> {log}
		hisat2-build --ss data/reference/splice-sites.gtf --exon data/reference/exon-sites.gtf {input.ref} {input.ref} 2>> {log}
		"""
	
	
rule HISAT2align:
	input:
		fwd="data/reads/{sample}_1.fastq.gz",
		rev="data/reads/{sample}_2.fastq.gz",
		idx="data/reference/ht.index"
	output:
		"data/alignments/{sample}.bam"
	log:
		align = "logs/hisat2/{sample}_align.log",
		sort = "logs/samtools_sort/{sample}.log"
	params:
		ref=lambda wildcards:config['ref']['genome'],
	threads:12
	shell:
		"""
		hisat2 -x {params.ref} -1 {input.fwd} -2 {input.rev} -p {threads} -q --dta --rg-id {wildcards.sample} --rg SM:{wildcards.sample} 2> {log.align} | 
		samtools view -bS - | samtools sort - -o {output} 2> {log.sort}
		"""

rule indexbams:
     input:
        "data/alignments/{sample}.bam"
     output:
        "data/alignments/{sample}.bam.bai"
     log:
        "logs/samtools_index/{sample}.log"
     shell:
        "samtools index {input} {output} 2> {log}"

rule mpileup:
    input:
        bam="data/alignments/{sample}.bam",
        idx="data/alignments/{sample}.bam.bai"
    output:
        "analysis/allele_balance/counts/{sample}_{mut}_allele_counts.tsv"
    log:
        "logs/mpileup/{sample}_{mut}.log"
    params:
        region = lambda wildcards: mutation_data[mutation_data.Name == wildcards.mut].Location.tolist(),
        ref = lambda wildcards:config['ref']['genome']
    shell:
        """
		mkdir -pv analysis/allele_balance/counts
        samtools mpileup {input.bam} -r {params.region} -f {params.ref} | python2 analysis/scripts/baseParser.py > {output}
        """

rule allele_balance:
    input:
        expand("analysis/allele_balance/counts/{sample}_{mut}_allele_counts.tsv", sample=samples, mut=mutations)
    output:
        "analysis/allele_balance/allele_balance.xlsx"
    log:
        "logs/allele_balance/Rscript.log"
    shell:
        """
		mkdir -pv analysis/allele_balance/csvs
		Rscript analysis/scripts/allele_balance.R 2> {log} 
        """


rule generateParams_Freebayes:
	input:
		ref_idx=lambda wildcards:config['ref']['genome'],
		metadata="data/samples.tsv"
	output:
		regions="data/regions/freebayes.regions",
		bamlist="data/bam.list",
		pops="data/populations.tsv"
	shell:
		"""
		ls data/alignments/*bam > {output.bamlist}
		cut -f 4,7 {input.metadata} | tail -n +2 > {output.pops}
		"""

rule variantCalling_Freebayes:
	input:
		bams=expand("data/alignments/{sample}.bam", sample=samples),
	    idx=expand("data/alignments/{sample}.bam.bai", sample=samples),
	    pops="data/populations.tsv",
		bamlist="data/bam.list"
	output:
		"data/variants/variants.freebayes.{chrom}.vcf"
	log:
		"logs/freebayes/{chrom}.log"
	params:
		ref=lambda wildcards:config['ref']['genome'],
		ploidy=10
	shell:
		"""
		freebayes -L {input.bamlist} -f {params.ref} \
		-r {wildcards.chrom} --populations {input.pops} --ploidy {params.ploidy} --pooled-discrete --use-best-n-alleles 5 > {output} 2> {log}
		"""


rule snpEff:
    input:
        vcf="data/variants/merged.vcf.gz",
	    database="Aedes_aegypti_lvpagwg"
    output:
        "data/variants/snpEff.merged.vcf.gz"
    log:
        "logs/SnpEff/snpeff.log"
    params:
        prefix="data/variants/snpEff.merged.vcf"
    shell:
        """
        java -jar ~/apps/snpEff/snpEff.jar download {input.database}
        java -jar ~/apps/snpEff/snpEff.jar eff {input.database} {input.vcf} > {params.prefix}
        bgzip {params.prefix}
        """