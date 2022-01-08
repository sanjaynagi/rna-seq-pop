################################          Common functions           ##################################

## If PBS is activated
if config['VariantAnalysis']['selection']["pbs"]["activate"]:
    windowedStats = ["Fst", "Pbs"]
else:
    windowedStats = ["Fst"]


def getFASTQs(wildcards, rules=None):
    """
    Get FASTQ files from unit sheet.
    If there are more than one wildcard (aka, sample), only return one fastq file
    If the rule is HISAT2align, then return the fastqs with -1 and -2 flags
    """
    
    if config["fastq"]["auto"]:
        units = pd.read_csv(config["metadata"], sep="\t")
        units = (
            units.assign(fq1=f"resources/reads/" + units["sampleID"] + "_1.fastq.gz")
            .assign(fq2=f"resources/reads/" + units["sampleID"] + "_2.fastq.gz")
            .set_index("sampleID")
        )
    else:
        assert os.path.isfile(
            config["fastq"]["table"]
        ), f"config['fastq']['table'] (the config/fastq.tsv file) does not seem to exist. Please create one, or use the 'auto' option and name the fastq files as specified in the config/README.md"
        units = pd.read_csv(config["fastq"]["table"], sep="\t", index_col="sampleID")

    if len(wildcards) > 1:
        u = units.loc[wildcards.sample, f"fq{wildcards.n}"]
        return u
    else:
        u = units.loc[wildcards.sample, ["fq1", "fq2"]].dropna()
        if rules == "HISAT2align":
            return [f"-1 {u.fq1} -2 {u.fq2}"]
        else:
            return [f"{u.fq1}", f"{u.fq2}"]


def GetDesiredOutputs(wildcards):

    """
    Function that returns a list of the desired outputs for the rule all, depending on the config.yaml
    configuration file. As of V0.4.0 Does not list every single output, but should mean all rules and desired outputs are created.
    """

    wanted_input = []

    # QC & Coverage
    if config["QualityControl"]["activate"]:
        wanted_input.extend(
            expand(
                [
                    "results/.input.check",
                    "resources/reads/qc/{sample}_{n}_fastqc.html",
                    "results/multiQC.html",
                ],
                sample=samples,
                n=[1, 2],
                chrom=config["chroms"],
            )
        )
        if config["VariantAnalysis"]["activate"]:
            wanted_input.extend(
            expand(
                [
                    "results/alignments/coverage/{sample}.mosdepth.summary.txt",
                    "results/alignments/bamStats/{sample}.flagstat",
                ],
                sample=samples,
            )
        )   
    
    if config['DifferentialExpression']['activate']:
        
        # Differential Expression outputs
        wanted_input.extend(
            expand(
                [
                    "results/genediff/{comp}.csv",
                    "results/genediff/{dataset}_diffexp.xlsx",
                    "results/isoformdiff/{comp}.csv",
                    "results/isoformdiff/{dataset}_isoformdiffexp.xlsx",
                    "results/plots/PCA.pdf",
                    "results/quant/countStatistics.tsv",
                ],
                comp=config["contrasts"],
                dataset=config["dataset"],
            )
        )

    if config['DifferentialExpression']["progressiveGenes"]["activate"]:
        wanted_input.extend(
            expand(
                "results/genediff/{name}.{direction}.progressive.tsv",
                name=config['DifferentialExpression']["progressiveGenes"]["groups"],
                direction=["up", "down"],
            )
        )

    if config["VariantAnalysis"]["activate"]:
        wanted_input.extend(
            expand(
                [
                    "results/variantAnalysis/vcfs/stats/{chrom}.txt",
                    "results/variantAnalysis/pca/PCA-{chrom}-{dataset}.png",
                    "results/variantAnalysis/SNPstats/snpsPerGenomicFeature.tsv",
                    "results/variantAnalysis/SNPstats/nSNPsPerGene.tsv",
                    "results/variantAnalysis/diversity/{dataset}_SNPdensity_{chrom}.png",
                    "results/variantAnalysis/diversity/SequenceDiversity.tsv",
                    "results/variantAnalysis/selection/FstPerGene.tsv",
                    "results/variantAnalysis/selection/TajimasDPerGene.tsv",
                    "results/variantAnalysis/diversity/SequenceDivPerGene.tsv",
                    "results/variantAnalysis/diversity/DxyPerGene.tsv",
                    "results/variantAnalysis/selection/fst/Fst.{comp}.{wsize}.{chrom}.png",
                ],
                chrom=config["chroms"],
                dataset=config["dataset"],
                comp=config["contrasts"],
                wsize=config['VariantAnalysis']['selection']["pbs"]["windownames"],
            )
        )


        if config['VariantAnalysis']['ploidy'] > 1:
            wanted_input.extend(
                expand(
                    [
                        "results/variantAnalysis/diversity/inbreedingCoef.tsv",
                    ]
                )
            )


    if config['VariantAnalysis']["AIMs"]["activate"]:
        wanted_input.extend(
            expand(
                [
                    "results/variantAnalysis/AIMs/AIMs_summary.tsv",
                    "results/variantAnalysis/AIMs/AIM_fraction_whole_genome.png",
                    "results/variantAnalysis/AIMs/n_AIMS_per_chrom.tsv",
                    "results/variantAnalysis/AIMs/AIM_fraction_{chrom}.tsv",
                ],
                chrom=config["chroms"],
            )
        )

    if config['miscellaneous']["VariantsOfInterest"]["activate"]:
        wanted_input.extend(
            [
                "results/variantAnalysis/variantsOfInterest/alleleBalance.xlsx",
                "results/variantAnalysis/variantsOfInterest/VOI.heatmapPerSample.png",
                "results/variantAnalysis/variantsOfInterest/VOI.heatmapPerTreatment.png"
                ]
            )

    if config['DifferentialExpression']["GSEA"]["activate"]:
        wanted_input.extend(
            expand(
                [
                    "results/gsea/genediff/{comp}.DE.{pathway}.tsv",
                ],
                comp=config["contrasts"],
                pathway=["kegg", "GO"],
            )
        )

    if config["VariantAnalysis"]["activate"] and config['DifferentialExpression']["GSEA"]["activate"]:
        wanted_input.extend(
            expand(
                ["results/gsea/fst/{comp}.FST.{pathway}.tsv"],
                comp=config["contrasts"],
                pathway=["kegg", "GO"],
            )
        )

    if config['miscellaneous']["diffsnps"]["activate"]:
        wanted_input.extend(
            expand(
                [
                    "results/variantAnalysis/diffsnps/{comp}.sig.kissDE.tsv",
                ],
                comp=config["contrasts"],
            )
        )

    if config['miscellaneous']["venn"]["activate"]:
        wanted_input.extend(
            expand(
                [
                    "results/RNA-Seq-full.xlsx",
                    "results/venn/{comp}_DE.Fst.venn.png",
                ],
                comp=config["contrasts"],
            )
        )

    if config['VariantAnalysis']["karyotype"]["activate"]:
        wanted_input.extend(
            expand(
                [
                    "results/karyotype/{karyo}.karyo.txt",
                    "results/karyotype/karyoFreqs.png"
                ],
                karyo=config['VariantAnalysis']["karyotype"]["inversions"],
            )
        )

    if config['miscellaneous']["sweeps"]["activate"]:
        if config['DifferentialExpression']['activate']:
            wanted_input.extend(
                expand(
                    [
                        "results/genediff/ag1000gSweeps/{comp}_swept.tsv",
                    ],
                    comp=config["contrasts"],
                )
            )

    # wanted_input.extend(["results/quant/percentageContributionGeneCategories.tsv"])

    return(wanted_input)


def welcome(version):

    print("---------------------------- RNA-Seq-Pop ----------------------------")
    print(f"Running RNA-Seq-Pop snakemake workflow in {workflow.basedir}\n")
    print(f"Author:   Sanjay Curtis Nagi")
    print(f"Workflow Version: {version}")
    print("Execution time: ", datetime.datetime.now())
    print(f"Dataset: {config['dataset']}", "\n")
