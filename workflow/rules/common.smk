import pandas as pd

def get_desired_outputs(wildcards): 

    """
    Function that returns a list of the desired outputs for the rule all, depending on the config.yaml
    configuration file. As of V0.1 Does not list every single output, but will mean all rules and outputs are created.
    """

    wanted_input = []

    # QC & Coverage
    if config['qc']['activate']:
        wanted_input.extend(
            expand(
            [
                "resources/reads/qc/{sample}_{n}_fastqc.html",
                "resources/alignments/coverage/{sample}.mosdepth.summary.txt",
            ],
            sample=samples, 
            n=[1,2]
            )
        )

    # Differential Expression outputs
    wanted_input.extend(
                [
                    "results/genediff/RNA-Seq_diff.xlsx",
                    "results/isoformdiff/RNA-Seq_isoformdiff.xlsx",
                    "results/plots/PCA.pdf",
                    "results/quant/count_statistics.tsv"
                ],
        )

    # Progressive changes in expression - e.g Kisumu < control field < resistant field
    if config['progressiveGenes']['activate']:
        wanted_input.extend(
            expand("results/genediff/{name}.{direction}.progressive.tsv", name=config['progressiveGenes']['groups'], direction=["up", "down"])
        )

    # Variant analyses
    wanted_input.extend(
        expand(
            [
                "results/variants/diffsnps/{name}.sig.kissDE.tsv",
                "results/variants/plots/PCA-{chrom}-{dataset}.png",
                "results/variants/plots/{dataset}_SNPdensity_{chrom}.png",
                "results/variants/stats/inbreedingCoef.tsv",
                "results/variants/stats/inbreedingCoef.mean.tsv",
                "results/variants/stats/SequenceDiversity.tsv",
                "results/variants/Fst.tsv",
                "results/variants/TajimasD.tsv",
                "results/variants/SequenceDiv.tsv",
                "results/variants/plots/fst/{comp}.{chrom}.fst.line.png",
        #        "results/variants/plots/pbs/{pbscomp}.{chrom}.pbs.line.png",
                "results/RNA-Seq-full.xlsx",
                "results/venn/{name}_DE.Fst.venn.png",
                ],
                name = config['contrasts'],
                chrom = config['chroms'],
                dataset = config['dataset'],
                stat = ['fst', 'pbs'],
                comp = config['contrasts'],
        #        pbscomp = config['pbs']['contrasts'],
                plot = ['line', 'scatter']
            )
        )

    # Ancestry Informative Markers
    if config['AIMs']['activate']:
        wanted_input.extend(
            expand(
                [
                    "results/variants/AIMs/AIMs_summary.tsv",
                    "results/variants/AIMs/AIM_fraction_whole_genome.png",
                    "results/variants/AIMs/n_AIMS_per_chrom.tsv",
                    "results/variants/AIMs/AIM_fraction_{chrom}.tsv",
                ],
                chrom=config['chroms'],
            )
        )

    #IRmutations
    if config['IRmutations']['activate']:
        wanted_input.extend(["results/allele_balance/allele_balance.xlsx"])
    

    if config['GSEA']['activate']:
       wanted_input.extend(
           expand(
                [
                    "results/gsea/genediff/{comp}.DE.csv",
                    "results/gsea/fst/{comp}.FST.csv",
                    "results/gsea/diffsnps/{comp}.diffsnps.csv",
                    "results/gsea/pbs/{pbscomp}.PBS.csv",
                ],
                comp = config['contrasts'],
                pbscomp = config['pbs']['contrasts'])
        )

    return(wanted_input)
















#def check_chroms(gffpath, chroms, str_remove):
#            """
#            Check that the chroms in the gff match ours, at least after removing the chosen string (str_remove)
#            """
#            
#            gff = allel.gff3_to_dataframe(gffpath)
# #           gff['seqid'] = gff['seqid'].str.strip(str_remove)
#                
#            assert np.isin(chroms, gff.seqid).all(), f"All provided chromosome names ({chroms}) are not present in GFF file"
#    

#rule check_inputs:
#    output:
##       touch("config/.input.check")
#    params:
# #      metadata=config['samples'],
#  #     chroms=config['chroms'],
#       gffpath=config['ref']['gff'],
#   ##    str_remove=config['ref']['str_remove'],
#       gene_names=config['ref']['genenames']
#    run:
#        metadata = pd.read_csv(params.metadata, sep="\t")
#        
#        #check to see if fastq files 
#        for sample in metadata.samples:
#            for n in [1,2]:
#                fqpath = f"resources/reads/{sample}_{n}.fastq.gz"
##
#                assert os.path.isfile(fqpath), "all sample names in 'samples.tsv' do not match a .fastq.gz file"#
#
# #       #check that the chromosomes are present in the gff file 
#  #      check_chroms(params.gffpath, params.chroms, params.str_remove)
#
# ##       #check column names of gene_names.tsv
#   #     gene_names = pd.read_csv(params.gene_names, sep="\t")
#    #    colnames = ['Gene_stable_ID', 'Gene_description', 'Gene_name']
#    #    assert (gene_names.columns == colnames).all(), f"Column names of gene_names.tsv should be {colnames}"


rule DElist:
    output:
        DElist="resources/DE.contrast.list"
    params:
        contrasts = config['contrasts']
    run:
        df = pd.DataFrame(params.contrasts)
        df.columns = ['contrast']
        df.to_csv(output.DElist, sep="\t", index=None)