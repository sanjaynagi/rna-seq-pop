################## lncrna ###############


rule stringtie:
    input:
        bam="data/alignments/{sample}.bam",
        gtf=lambda wildcards: config["ref"]["gtf"],
    output:
        "data/transcript_assembly/{sample}.gtf",
    log:
        "logs/stringtie/{sample}.log",
    threads: 12
    shell:
        """
        stringtie {input.bam} -G {input.gtf} -p {threads} -o {output} 2> {log}
        """


rule stringtie_merge:
    input:
        gtfs=expand("data/transcript_assembly/{sample}.gtf", sample=samples),
        refgtf=lambda wildcards: config["ref"]["gtf"],
    output:
        merged="data/transcript_assembly/merged.gtf",
        gtflist="data/transcript_assembly/all_gtfs.txt",
    log:
        "logs/stringtie/merge.log",
    threads: 12
    shell:
        """
        find data/transcript_assembly/ -maxdepth 1 -type f > all_gtfs.txt 2> {log}
        stringtie --merge -p {threads} -G {input.refgtf} -o {output.merged} all_gtfs.txt 2>> {log}
                mv all_gtfs.txt {output.gtflist}
        """


# filter gtf to > 200 bases and multi-exon
# gffread merged.annotated.gtf -l 200 -U > filtered.merged.gtf


### get transcripts fasta
# gffread -w transcripts.fa -g ../reference/aedes-aegypti-lvpagwg_chromosomes.l5.fa filtered.merged.gtf

# plek
# PLEK.py -fasta transcripts.fa -out transcripts_predicted -thread 10 -isoutmsg 1
