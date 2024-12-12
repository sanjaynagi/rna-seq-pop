
rule HISAT2index:
    """
    Make a HISAT2 index of the reference genome
    """
    input:
        fasta=config["reference"]["genome"].rstrip(".gz"),
    output:
        "resources/reference/ht2index/idx.1.ht2",
        touch("resources/reference/ht2index/.complete"),
    log:
        "logs/HISAT2/HISAT2index.log",
    conda:
        "../envs/variants.yaml"
    params:
        prefix=lambda w, output: output[0].split(os.extsep)[0],
    threads: 8
    shell:
        "hisat2-build -p {threads} {input.fasta} {params.prefix}  2> {log}"


rule HISAT2align:
    """
    Align reads to the genome with HISAT2, mark duplicates with samblaster and sort with samtools
    """
    input:
        reads=lambda wildcards: getFASTQs(
            wildcards=wildcards, rules="HISAT2align_input"
        ),
        idx="resources/reference/ht2index/.complete",
    output:
        "results/alignments/{sample}.hisat2.bam",
    log:
        align="logs/HISAT2/{sample}_align.log",
        sort="logs/samtoolsSort/{sample}.log",
    conda:
        "../envs/variants.yaml"
    params:
        readflags=lambda wildcards: getFASTQs(wildcards=wildcards, rules="HISAT2align"),
        extra="--dta -q --rg-id {sample} --rg SM:{sample} --rg PL:ILLUMINA --new-summary",
        idx="resources/reference/ht2index/idx",
        samblaster="" if config['fastq']['paired'] is True else "--ignoreUnmated"
    threads: 12
    shell:
        """
        hisat2 {params.extra} --threads {threads} -x {params.idx} {params.readflags} 2> {log.align} | 
        samblaster {params.samblaster} 2> {log.sort} | samtools sort -@{threads} - -o {output} 2>> {log.sort}
        """


chunks = np.arange(1, config["VariantAnalysis"]["chunks"])

rule GenerateFreebayesParams:
    input:
        ref_idx=config["reference"]["genome"].rstrip(".gz"),
        index=config["reference"]["genome"].rstrip(".gz") + ".fai",
        bams=expand("results/alignments/{sample}.hisat2.bam", sample=samples),
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
        metadata['bams'] = 'results/alignments/' + metadata['sampleID'] + '.hisat2.bam'
        
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


rule VariantCallingFreebayes:
    """
    Run freebayes on chunks of the genome, splitting the samples by population (strain)
    """
    input:
        bams=expand("results/alignments/{sample}.hisat2.bam", sample=samples),
        index=expand("results/alignments/{sample}.hisat2.bam.bai", sample=samples),
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


