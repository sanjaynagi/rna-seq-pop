################ 	Variant analysis 	##################
# Scikit-allel () ,                                               		
# samtools mpileup (Li et al., 2009)                                             
# kissDE (Kim et al., 2018)   

rule mpileupIR:
    input:
        bam="resources/alignments/{sample}.bam",
        index="resources/alignments/{sample}.bam.bai"
    output:
        "results/allele_balance/counts/{sample}_{mut}_allele_counts.tsv"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/mpileup/{sample}_{mut}.log"
    params:
        region = lambda wildcards: mutation_data[mutation_data.Name == wildcards.mut].Location.tolist(),
        ref = config['ref']['genome']
    shell:
        """
        samtools mpileup {input.bam} -r {params.region} -f {params.ref} | python2 workflow/scripts/baseParser.py > {output}
        """

rule alleleBalanceIR:
    input:
        counts = expand("results/allele_balance/counts/{sample}_{mut}_allele_counts.tsv", sample=samples, mut=mutations),
        samples = config['samples'],
        mutations = config['IRmutations']['path']
    output:
        allele_balance = "results/allele_balance/allele_balance.xlsx",
        mean_allele_balance = "results/allele_balance/mean_allele_balance.xlsx"
    conda:
        "../envs/alleleBalance.yaml"
    log:
        "logs/AlleleBalance.log"
    script:
        "../scripts/AlleleBalance.R"

rule alleleTables:
    input:
        bam="resources/alignments/{sample}.bam",
        bed="resources/regions/missense.pos.{chrom}.bed",
        ref=config['ref']['genome']
    output:
        "results/variants/alleleTables/{sample}.chr{chrom}.allele.table"
    conda:
        "../envs/variants.yaml"
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

rule DifferentialSNPs:
    input:
        samples=config['samples'],
        gff=config['ref']['gff'],
        DEcontrasts="resources/DE.contrast.list",
        geneNames = "resources/gene_names.tsv",
        tables=expand("results/variants/alleleTables/{sample}.chr{chrom}.allele.table", sample=samples, chrom=config['chroms'])
    output:
        expand("results/variants/diffsnps/{name}.sig.kissDE.tsv", name = config['contrasts']),
        expand("results/variants/diffsnps/{name}.kissDE.tsv", name = config['contrasts']),
        expand("results/variants/diffsnps/{name}.normcounts.tsv", name = config['contrasts'])
    conda:
        "../envs/diffsnps.yaml"
    log:
        "logs/variants/kissDE.log"
    params:
        chroms = config['chroms'],
        gffchromprefix="AaegL5_", # in case like the aedes genome, there is an annoying prefix before each chromosome
        mincounts = 100,
        pval_flt = 0.001, # pvalues already adjusted but way want extra filter for sig file
        prefix = "AaegL5_"
    script:
        "../scripts/differentialSNPs.R"

rule pca_variantdensity:
    input:
        vcf=expand("results/variants/annot.variants.{chrom}.vcf.gz", chrom=config['chroms'])
    output:
        pcafig=expand("results/variants/plots/PCA-{chrom}-{dataset}.png", chrom=config['chroms'], dataset=config['dataset']),
        snpdensityfig=expand("results/variants/plots/{dataset}_SNPdensity_{chrom}.png", chrom=config['chroms'], dataset=config['dataset'])
    log:
        "logs/pca/pca.log"
    conda:
        "../envs/fstpca.yaml"
    params:
        dataset = config['dataset'],
        chroms = config['chroms'],
        missingprop = 0.98
    script:
        "../scripts/pca.py"

rule Fst_PBS_TajimaD_SeqDiv_per_gene:
    input:
        samples=config['samples'],
        gff=config['ref']['gff'],
        DEcontrasts="resources/DE.contrast.list",
        geneNames = "resources/gene_names.tsv",
        vcf=expand("results/variants/annot.variants.{chrom}.vcf.gz", chrom=config['chroms'])
    output:
        "results/variants/Fst_PBS.tsv",
        "results/variants/tajimas_d.tsv",
        "results/variants/sequence_div.tsv"
    conda:
        "../envs/fstpca.yaml"
    log:
        "logs/variants/fst_pbs.log"
    params:
        pbs = config['pbs']['activate'],
        pbscomps = config['pbs']['contrasts'],
        chroms = config['chroms'],
        ploidy = config['ploidy'],
        missingprop = 0.66,
        gffchromprefix=config['ref']['str_remove'] #in case like the aedes genome, there is an annoying before each chromosome
    script:
        "../scripts/FstPbsTajimasDseqDiv.py"

rule Venn:
   input:
        DEcontrasts="resources/DE.contrast.list",
        DE="results/genediff/RNA-Seq_diff.xlsx",
        FST="results/variants/Fst_PBS.tsv",
        diffsnps=expand("results/variants/diffsnps/{name}.sig.kissDE.tsv", name = config['contrasts'])
   output:
        "results/RNA-Seq-full.xlsx",
        expand("results/venn/{name}_DE.Fst.venn.png", name=config['contrasts'])
   conda:
        "../envs/fstpca.yaml"
   log:
        "logs/venn.log"
   params:
        pbs=config['pbs']['activate'],
        pbscomps=config['pbs']['contrasts'],
        percentile=0.05
   script:
       "../scripts/venn.py"
