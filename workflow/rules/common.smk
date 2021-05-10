################################          Common functions           ##################################

## If PBS is activated 
if config['pbs']:
    windowedStats = ['fst', 'pbs']
else:
    windowedStats = ['fst']


def getFASTQs(wildcards):
    """
    Get FASTQ files from unit sheet.
    If there are more than one wildcard (aka, sample), only return one fastq file
    """
    if config['fastq']['auto']:
        units = pd.read_csv(config['samples'], sep="\t")
        units = units[['samples']].assign(fq1=units['samples'] + "_1.fastq.gz").assign(fq2=units['samples'] + "_2.fastq.gz").set_index('samples')
    else:
        assert os.path.isfile(config['fastq']['table']), f"config['fastq']['table'] (the config/fastq.tsv file) does not seem to exist. Please create one, or use the 'auto' option and name the fastq files as specified in the config/README.md"
        units = pd.read_csv(config['fastq']['table'], sep="\t", index_col="samples")

    if len(wildcards) > 1:
        u = units.loc[wildcards.sample, f"fq{wildcards.n}"]
        return u
    else:
        u = units.loc[wildcards.sample, ["fq1", "fq2"]].dropna()
        return [f"{u.fq1}", f"{u.fq2}"]


def get_desired_outputs(wildcards): 

    """
    Function that returns a list of the desired outputs for the rule all, depending on the config.yaml
    configuration file. As of V0.4.0 Does not list every single output, but should mean all rules and desired outputs are created.
    """

    wanted_input = []

    # QC & Coverage
    if config['qc']['activate']:
        wanted_input.extend(
            expand(
            [
                "config/.input.check",
                "resources/reads/qc/{sample}_{n}_fastqc.html",
                "resources/alignments/coverage/{sample}.mosdepth.summary.txt",
                "resources/alignments/bamStats/{sample}.flagstat",
                "results/multiQC.html"
            ],
            sample=samples, 
            n=[1,2],
            chrom = config['chroms']
            )
        )

    # Differential Expression outputs
    wanted_input.extend(
            expand(
                [
                    "results/genediff/{comp}.csv",
		            "results/genediff/{dataset}_diffexp.xlsx",
                    "results/isoformdiff/{comp}.csv",
                    "results/isoformdiff/{dataset}_isoformdiffexp.xlsx",
                    "results/plots/PCA.pdf",
                    "results/quant/count_statistics.tsv"
                ],
             comp = config['contrasts'],
             dataset = config['dataset']
            )
        )

    if config['progressiveGenes']['activate']:
        wanted_input.extend(
            expand("results/genediff/{name}.{direction}.progressive.tsv", name=config['progressiveGenes']['groups'], direction=["up", "down"])
        )

    if config['VariantCalling']['activate']:
        wanted_input.extend(
            expand(
                [   
                    "results/variants/vcfs/stats/{chrom}.txt",
                    "results/variants/plots/PCA-{chrom}-{dataset}.png",
                    "results/variants/plots/{dataset}_SNPdensity_{chrom}.png",
                    "results/variants/stats/inbreedingCoef.tsv",
                    "results/variants/stats/inbreedingCoef.mean.tsv",
                    "results/variants/stats/SequenceDiversity.tsv",
                    "results/variants/fst.tsv",
                    "results/variants/TajimasD.tsv",
                    "results/variants/SequenceDiv.tsv",
                    "results/variants/plots/fst/{comp}.{chrom}.fst.{wsize}.png",
                    ],
                    chrom = config['chroms'],
                    dataset = config['dataset'],
                    comp = config['contrasts'],
                    wsize = config['pbs']['windownames']
                )
            )

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

    if config['IRmutations']['activate']:
        wanted_input.extend(["results/allele_balance/allele_balance.xlsx"])
    

    if config['GSEA']['activate']:
       wanted_input.extend(
           expand(
               [
               "results/gsea/genediff/{comp}.DE.{pathway}.tsv",
                ],
                comp = config['contrasts'],
                pathway=['kegg', 'GO']
                )
        )
    
    if config['VariantCalling']['activate'] and config['GSEA']['activate']:
        wanted_input.extend(
            expand(
            [
                "results/gsea/fst/{comp}.FST.{pathway}.tsv"
            ],
            comp = config['contrasts'],
            pathway=['kegg', 'GO']
            )
        )


    if config['diffsnps']['activate']:
       wanted_input.extend(
           expand(
               [
                "results/variants/diffsnps/{comp}.sig.kissDE.tsv",
                ],
                comp = config['contrasts'],
                )
	    )

    if config['venn']['activate']:
       wanted_input.extend(
           expand(
               [
                "results/RNA-Seq-full.xlsx",
                "results/venn/{comp}_DE.Fst.venn.png",
                ],
                comp = config['contrasts'],
                )
	    )

    if config['karyotype']['activate']:
        wanted_input.extend(
            expand(
                [
                "results/karyotype/{karyo}.karyo.txt"
                ],
                karyo=config['karyotype']['inversions']
            )
        )
    
    if config['sweeps']['activate']:
        wanted_input.extend(
            expand(
                [
                    "results/genediff/ag1000gSweeps/{comp}_swept.tsv",
                ],
                comp=config['contrasts']
            )
        )

    #wanted_input.extend(["results/quant/percentageContributionGeneCategories.tsv"])

    return(wanted_input)


def welcome(version, verbose=False):

    print("---------------------------- RNA-Seq-IR ----------------------------")
    print(f"Running RNA-Seq-IR snakemake workflow in {workflow.basedir}\n")
    print(f"Author:   Sanjay Curtis Nagi")
    print(f"Workflow Version: {version}")
    print("Execution time: ", datetime.datetime.now())
    print(f"Dataset: {config['dataset']}", "\n")
    
    if verbose:
        print("Modules activated:", "\n")
            
        if config['qc']['activate']: 
            print('QC - FASTQC, mosdepth, samtools flagstat')

        print("Differential Expression - Kallisto, DESeq2, Sleuth")
        print("Writing results to results/genediff, results/isoformdiff")
        print(f"contrasts: {config['contrasts']}\n")

        if config['progressiveGenes']['activate']:
            print("ProgressiveGenes - Checking for DE genes differentially expressed in the same direction")
            print(f"groups: {config['progressiveGenes']['groups']}\n")

        if config['GSEA']['activate']:
            print("GSEA")

        if config['IRmutations']['activate']:
            print(f"IRmutations - reporting raw allele frequencies from loci found in {config['IRmutations']['path']}")
        
        if config['pbs']['activate']:
            print("PBS - Performing population branch statistic analysis")
            print(f"contrasts: {config['pbs']['contrasts']}", "\n")
        

        if config['diffsnps']['activate']:
            print("DifferentialSNPs - Performing differential SNP testing analysis")
        
        if config['AIMs']['activate']:
            print("AIMs - Ancestry informative marker analysis")