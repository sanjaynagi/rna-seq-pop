chunks = np.arange(1, config["VariantAnalysis"]["chunks"])

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


rule RestrictToSNPs:
    """"
    Filter out indels
    """
    input:   
        "results/variantAnalysis/vcfs/freebayes/variants.{contig}.vcf" if config['VariantAnalysis']['caller'] == 'freebayes' else "results/variantAnalysis/vcfs/haplotypecaller/variants.{contig}.vcf",
    output:
        temp("results/variantAnalysis/vcfs/variants.{contig}.vcf"),
    log:
        "logs/bcftoolsView/{contig}.log",
    params:
        extra="-v snps",
    wrapper:
        "v1.15.0/bio/bcftools/view"


rule snpEffDbDownload:
    """
    Download the snpEff database for your species
    """
    output:
        touch("resources/reference/mysnpeffdb/db.dl"),
    log:
        "logs/snpEff/snpEffDbDownload.log",
    conda:
        "../envs/snpeff.yaml"
    params:
        ref=config["reference"]["snpeff"]['dbname'],
        dir="resources/reference/mysnpeffdb",
    shell:
        "snpEff download {params.ref} -dataDir {params.dir} 2> {log}"

rule createCustomSnpEffDb:
    """
    Create a custom SnpEff database from a reference genome and GFF file
    """
    input:
        fa="path/to/your_organism.fa",
        gff="path/to/your_organism.gff",
    output:
        touch("workflow/scripts/snpEff/custom_db.built"),
    log:
        "logs/snpEff/createCustomSnpEffDb.log",
    conda:
        "../envs/snpeff.yaml"
    params:
        db_name="my_organism",
        data_dir="resources/snpEffdb/my_organism",
    shell:
        """
        mkdir -p {params.data_dir}
        cp {input.fa} {params.data_dir}/sequences.fa
        cp {input.gff} {params.data_dir}/genes.gff
        echo "{params.db_name}.genome : My Organism" >> snpEff.config
        snpEff build -gff3 -v {params.db_name} 2> {log}
        touch {output}
        """


rule createCustomSnpEffDb:
    """
    Create a custom SnpEff database from a reference genome and GFF file
    """
    input:
        fa=config['reference']['genome'].rstrip(".gz"),
        gff=config['reference']['gff'],
    output:
        "resources/reference/mysnpeffdb/sequences.fa",
        "resources/reference/mysnpeffdb/genes.gff",
        "resources/reference/mysnpeffdb/snpEffectPredictor.bin"
    log:
        "logs/snpEff/createCustomSnpEffDb.log",
    conda:
        "../envs/snpeff.yaml"
    params:
        dataDir=lambda x: wkdir + "/resources/reference",
        wkdir=wkdir
    shell:
        """
        ln -s {params.wkdir}/{input.fa} {params.dataDir}/mysnpeffdb/sequences.fa 2> {log}
        ln -s {params.wkdir}/{input.gff} {params.dataDir}/mysnpeffdb/genes.gff 2> {log}
        snpEff build -gff3 -v -dataDir {params.dataDir} -configOption mysnpeffdb.genome=mysnpeffdb mysnpeffdb -noCheckCds -noCheckProtein 2>> {log}
        """


rule snpEff:
    """
    Run snpEff on the VCFs 
    """
    input:
        calls="results/variantAnalysis/vcfs/variants.{contig}.vcf",
        db="resources/reference/mysnpeffdb/snpEffectPredictor.bin" if config['reference']['snpeff']['customdb'] else "resources/reference/mysnpeffdb/db.dl",
    output:
        calls=expand(
            "results/variantAnalysis/vcfs/{dataset}.{{contig}}.vcf.gz",
            dataset=config["dataset"],
        ),
        csvStats="results/variantAnalysis/vcfs/snpEff.summary.{contig}.csv",
    log:
        "logs/snpEff/snpEff.{contig}.log",
    conda:
        "../envs/snpeff.yaml"
    params:
        db=config["reference"]["snpeff"]['dbname'],
        prefix=lambda w, output: os.path.splitext(output[0])[0],
        dataDir=lambda x: wkdir + "resources/reference/mysnpeffdb",
    shell:
        """
        snpEff eff {params.db} -dataDir {params.dataDir} -csvStats {output.csvStats} {input.calls} > {params.prefix} 2> {log}
        bgzip {params.prefix}
        """

