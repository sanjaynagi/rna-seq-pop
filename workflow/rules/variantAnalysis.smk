################ Variant Analysis ##################
# samtools
# scikit-allel (Miles, Harding) 10.5281/zenodo.3935797
# kiss-DE

rule mpileupIR:
    input:
        bam = "resources/alignments/{sample}.bam",
        index = "resources/alignments/{sample}.bam.bai"
    output:
        "results/allele_balance/counts/{sample}_{mut}_allele_counts.tsv"
    conda:
        "../envs/variants.yaml"
    priority: 10
    log:
        "logs/mpileupIR/{sample}_{mut}.log"
    params:
        region = lambda wildcards: mutationData[mutationData.Name == wildcards.mut].Location.tolist(),
        ref = config['ref']['genome']
    shell:
        """
        samtools mpileup {input.bam} -r {params.region} -f {params.ref} | python2 {workflow.basedir}/scripts/BaseParser.py > {output}
        """

rule AlleleBalanceIR:
    input:
        counts = expand("results/allele_balance/counts/{sample}_{mut}_allele_counts.tsv", sample=samples, mut=mutationData.Name),
        samples = config['samples'],
        mutations = config['IRmutations']['path']
    output:
        expand("results/allele_balance/csvs/{mut}_allele_balance.csv", mut=mutationData.Name),
        allele_balance = "results/allele_balance/allele_balance.xlsx",
        mean_allele_balance = "results/allele_balance/mean_allele_balance.xlsx"
    conda:
        "../envs/diffexp.yaml"
    priority: 10
    log:
        "logs/AlleleBalanceIR.log"
    script:
        "../scripts/MutAlleleBalance.R"

rule AlleleTables:
    input:
        bam = "resources/alignments/{sample}.bam",
        bed = "resources/regions/missense.pos.{chrom}.bed",
        ref = config['ref']['genome']
    output:
        "results/variants/alleleTables/{sample}.chr{chrom}.allele.table"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/AlleleTables/{sample}.{chrom}.log"
    params:
        ignore_indels='false',
        baseflt=-5,
        min_alt=3,
        min_af=0
    shell:
        """
        samtools mpileup -f {input.ref} -l {input.bed} {input.bam} | 
        {workflow.basedir}/scripts/mpileup2readcounts/mpileup2readcounts 0 {params.baseflt} {params.ignore_indels} {params.min_alt} {params.min_af} > {output} 2> {log}
        """

rule DifferentialSNPs:
    input:
        samples = config['samples'],
        gff = config['ref']['gff'],
        DEcontrasts = "resources/DE.contrast.list",
        geneNames = "resources/gene_names.tsv",
        tables = expand("results/variants/alleleTables/{sample}.chr{chrom}.allele.table", sample=samples, chrom=config['chroms'])
    output:
        expand("results/variants/diffsnps/{name}.sig.kissDE.tsv", name = config['contrasts']),
        expand("results/variants/diffsnps/{name}.kissDE.tsv", name = config['contrasts']),
        expand("results/variants/diffsnps/{name}.normcounts.tsv", name = config['contrasts'])
    conda:
        "../envs/diffsnps.yaml"
    log:
        "logs/DifferentialSNPs.log"
    params:
        chroms = config['chroms'],
        mincounts = 100,
        pval_flt = 0.001, # pvalues already adjusted but way want extra filter for sig file
    script:
         "../scripts/DifferentialSNPs.R"


rule WindowedStatisticsAndPCA:
    input:
        vcf = expand("results/variants/vcfs/annot.variants.{chrom}.vcf.gz", chrom=config['chroms']),
        samples = config['samples'],
        contrasts = "resources/DE.contrast.list",
        gff = config['ref']['gff']
    output:
        PCAfig = expand("results/variants/plots/PCA-{chrom}-{dataset}.png", chrom=config['chroms'], dataset=config['dataset']),
        SNPdensityFig = expand("results/variants/plots/{dataset}_SNPdensity_{chrom}.png", chrom=config['chroms'], dataset=config['dataset']),
        Fst = expand("results/variants/plots/fst/{comp}.{chrom}.fst.{wsize}.png", comp=config['contrasts'], chrom=config['chroms'], wsize=config['pbs']['windownames']),
        PBS = expand("results/variants/plots/pbs/{pbscomp}.{chrom}.pbs.{wsize}.png", pbscomp=config['pbs']['contrasts'], chrom=config['chroms'], wsize=config['pbs']['windownames']) if config['pbs'] else [], 
        inbreedingCoef = "results/variants/stats/inbreedingCoef.tsv",
        inbreedingCoefMean = "results/variants/stats/inbreedingCoef.mean.tsv",
        SequenceDiversity = "results/variants/stats/SequenceDiversity.tsv",
        #LD = "results/variants/stats/LD.tsv"
    log:
        "logs/WindowedStatisticsAndPCA.log"
    conda:
        "../envs/fstpca.yaml"
    params:
        dataset = config['dataset'],
        chroms = config['chroms'],
        ploidy = config['ploidy'],
        pbs = config['pbs']['activate'],
        pbscomps = config['pbs']['contrasts'],
        missingprop = config['pbs']['missingness'],
        qualflt = 30,
        linkage = False,
        windowsizes = config['pbs']['windowsizes'],
        windowsteps = config['pbs']['windowsteps'],
        windownames = config['pbs']['windownames']
    script:
        "../scripts/WindowedStatsAndPCA.py"

rule FstPbsPerGene:
    input:
        samples = config['samples'],
        gff = config['ref']['gff'],
        DEcontrasts = "resources/DE.contrast.list",
        geneNames = "resources/gene_names.tsv",
        vcf = expand("results/variants/vcfs/annot.variants.{chrom}.vcf.gz", chrom=config['chroms'])
    output:
        expand("results/variants/{stat}.tsv", stat=windowedStats),
        "results/variants/TajimasD.tsv",
        "results/variants/SequenceDiv.tsv"
    conda:
        "../envs/fstpca.yaml"
    log:
        "logs/FstPbsPerGene.log"
    params:
        pbs = config['pbs']['activate'],
        pbscomps = config['pbs']['contrasts'],
        chroms = config['chroms'],
        ploidy = config['ploidy'],
        missingprop = 0.8,
    script:
        "../scripts/FstPbsPerGene.py"


rule AncestryInformativeMarkers:
    input:
        vcf = expand("results/variants/vcfs/annot.variants.{chrom}.vcf.gz", chrom=config['chroms']),
        samples = config['samples'],
        aims_zarr_gambcolu = config['AIMs']['gambcolu'],
        aims_zarr_arab = config['AIMs']['arab']
    output:
        AIMs = "results/variants/AIMs/AIMs_summary.tsv",
        AIMs_fig = "results/variants/AIMs/AIM_fraction_whole_genome.png",
        n_AIMs = "results/variants/AIMs/n_AIMS_per_chrom.tsv",
        AIMs_chroms = expand("results/variants/AIMs/AIM_fraction_{chrom}.tsv", chrom=config['chroms'])
    log:
        "logs/AncestryInformativeMarkers.log"
    conda:
        "../envs/fstpca.yaml"
    params:
        chroms = config['chroms'],
        ploidy = config['ploidy'],
        missingprop = config['AIMs']['missingness'],
        qualflt = 30
    script:
        "../scripts/AncestryInformativeMarkers.py"


rule Karyotype:
    input:
        vcf = lambda wildcards: "results/variants/vcfs/annot.variants.2L.vcf.gz" if wildcards.karyo == "2La" else "results/variants/vcfs/annot.variants.2R.vcf.gz"
    output:
        "results/karyotype/{karyo}.karyo.txt"
    log:
        "logs/compKaryo/Karyo_{karyo}.log"
    conda:
        "../envs/fstpca.yaml"
    params:
         ploidy = config['ploidy']
    shell:
        """
        paste <(bcftools query -l {input.vcf}) \
        <(python {workflow.basedir}/scripts/compkaryo/compkaryo/compkaryo.py {input.vcf} {wildcards.karyo} -p {params.ploidy}) | 
        column -s $'\\t' -t > {output}
        """

        
rule VennDiagrams:
   input:
        DEcontrasts = "resources/DE.contrast.list",
        DE = "results/genediff/RNA-Seq_diff.xlsx",
        Fst = expand("results/variants/{stat}.tsv", stat=windowedStats),
        diffsnps = expand("results/variants/diffsnps/{name}.sig.kissDE.tsv", name = config['contrasts']) if config['diffsnps']['activate'] else []
   output:
        "results/RNA-Seq-full.xlsx",
        expand("results/venn/{name}_DE.Fst.venn.png", name=config['contrasts'])
   conda:
        "../envs/fstpca.yaml"
   log:
        "logs/VennDiagrams.log"
   params:
        pbs = config['pbs']['activate'],
        pbscomps = config['pbs']['contrasts'],
        percentile = 0.05
   script:
       "../scripts/VennDiagrams.py"