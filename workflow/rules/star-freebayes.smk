
rule star_index_cpu:
    """
    Build a STAR genome index for CPU alignment.
    """
    input:
        fasta=config["reference"]["genome"].rstrip(".gz"),
        gff=config["reference"]["gff"],
    output:
        index=directory("resources/reference/star_index_cpu"),
    log:
        "logs/star/index_cpu.log",
    conda:
        "../envs/variants.yaml"
    params:
        extra="--sjdbGTFtagExonParentTranscript Parent",
    threads: 8
    shell:
        """
        mkdir -p {output.index}
        STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gff} \
            {params.extra} \
            --genomeDir {output.index} \
            --outFileNamePrefix {output.index}/ > {log} 2>&1
        """


rule star_align_cpu:
    """
    Align reads with STAR (CPU) and output a coordinate-sorted BAM plus index.
    """
    input:
        reads=lambda wildcards: get_fastqs(wildcards=wildcards, rules="star_align_cpu_input"),
        idx="resources/reference/star_index_cpu",
    output:
        bam="results/alignments/{sample}.star.bam",
        bai="results/alignments/{sample}.star.bam.bai",
    log:
        align="logs/star/{sample}_align_cpu.log",
        sort="logs/samtoolsSort/{sample}.log",
    conda:
        "../envs/variants.yaml"
    params:
        prefix="results/alignments/{sample}.star.",
    threads: 12
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.idx} \
            --readFilesIn {input.reads} \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within \
            --outSAMattributes Standard \
            --outFileNamePrefix {params.prefix} > {log.align} 2>&1
        mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam}
        samtools index -@ {threads} {output.bam} {output.bai} 2> {log.sort}
        """


chunks = np.arange(1, config["VariantAnalysis"]["chunks"])

rule generate_freebayes_params:
    input:
        ref_idx=config["reference"]["genome"].rstrip(".gz"),
        index=config["reference"]["genome"].rstrip(".gz") + ".fai",
        bams=expand("results/alignments/{sample}.star.bam", sample=samples),
    output:
        bamlist="results/alignments/bam.list",
        pops="results/alignments/populations.tsv",
        regions=expand(
            "results/variantAnalysis/regions/genome.{contig}.region.{i}.bed",
            contig=config["contigs"],
            i=chunks,
        ),
    log:
        "logs/GenerateFreebayesParams.log",
    params:
        metadata=config["metadata"],
        contigs=config["contigs"],
        chunks=config["VariantAnalysis"]["chunks"],
    run:
        import pandas as pd
        import numpy as np
        from pathlib import Path
        import logging
        
        # Setup logging
        logging.basicConfig(
            filename=log[0],
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        
        # Load and filter fai file
        fai = pd.read_csv(input.index, sep='\t', header=None, usecols=[0, 1])
        fai = fai[fai[0].isin(params.contigs)]
        
        # Create bed files for each contig and chunk
        for contig in params.contigs:
            contig_length = fai[fai[0] == contig][1].iloc[0]
            bedseq = np.round(np.linspace(0, contig_length, params.chunks))
            
            for i in range(params.chunks - 1):
                bed_content = f"{contig}\t{int(bedseq[i])}\t{int(bedseq[i+1])}"
                output_path = Path(f"results/variantAnalysis/regions/genome.{contig}.region.{i+1}.bed")
                output_path.parent.mkdir(parents=True, exist_ok=True)
                
                with open(output_path, 'w') as f:
                    f.write(bed_content)
        
        # Load metadata
        metadata_path = params.metadata
        file_extension = Path(metadata_path).suffix.lower()
        
        if file_extension == '.xlsx':
            metadata = pd.read_excel(metadata_path)
        elif file_extension == '.tsv':
            metadata = pd.read_csv(metadata_path, sep='\t')
        elif file_extension == '.csv':
            metadata = pd.read_csv(metadata_path)
        else:
            raise ValueError("Metadata file must be .xlsx, .tsv, or .csv")
        
        # Add bam paths and create output files
        metadata['bams'] = 'results/alignments/' + metadata['sampleID'] + '.star.bam'
        
        # Create populations file
        metadata[['bams', 'strain']].to_csv(
            output.pops, 
            sep='\t', 
            header=False, 
            index=False
        )
        
        # Create bamlist file
        metadata[['bams']].to_csv(
            output.bamlist, 
            sep='\t', 
            header=False, 
            index=False
        )


rule variant_calling_freebayes:
    """
    Run freebayes on chunks of the genome, splitting the samples by population (strain)
    """
    input:
        bams=expand("results/alignments/{sample}.star.bam", sample=samples),
        index=expand("results/alignments/{sample}.star.bam.bai", sample=samples),
        ref=config["reference"]["genome"].rstrip(".gz"),
        samples="results/alignments/bam.list",
        pops="results/alignments/populations.tsv",
        regions="results/variantAnalysis/regions/genome.{contig}.region.{i}.bed",
    output:
        temp("results/variantAnalysis/vcfs/freebayes/{contig}/variants.{i}.vcf"),
    log:
        "logs/VariantCallingFreebayes/{contig}.{i}.log",
    params:
        ploidy=config["VariantAnalysis"]["ploidy"],
    conda:
        "../envs/variants.yaml"
    threads: 1
    shell:
        "freebayes -f {input.ref} -t {input.regions} --ploidy {params.ploidy} --populations {input.pops} --pooled-discrete --use-best-n-alleles 5 --min-alternate-fraction 0.05 -L {input.samples} > {output} 2> {log}"

chunks = np.arange(1, config["VariantAnalysis"]["chunks"])

rule concat_vcfs:
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
