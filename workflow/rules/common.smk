if config["VariantAnalysis"]["selection"]["pbs"]["activate"]:
    windowedStats = ["Fst", "Pbs"]
else:
    windowedStats = ["Fst"]


def load_metadata(metadata_path):
    # load panel metadata
    if metadata_path.endswith('.xlsx'):
        metadata = pd.read_excel(metadata_path, engine='openpyxl')
    elif metadata_path.endswith('.tsv'):
        metadata = pd.read_csv(metadata_path, sep="\t")
    elif metadata_path.endswith('.csv'):
        metadata = pd.read_csv(metadata_path, sep=",")
    else:
        raise ValueError("Metadata file must be .xlsx or .csv")
    return metadata

    
rule set_kernel:
    input:
        f'{workflow.basedir}/envs/pythonGenomics.yaml'
    output:
        touch("results/.kernel.set")
    conda: f'{workflow.basedir}/envs/pythonGenomics.yaml'
    log:
        "logs/notebooks/set_kernel.log"
    shell: 
        """
        python -m ipykernel install --user --name=pythonGenomics 2> {log}
        """

rule remove_input_param_cell_nbs:
    input:
        input_nb = f"{workflow.basedir}/notebooks/process-notebooks.ipynb",
        counts_qc = "docs/rna-seq-pop-results/notebooks/counts-qc.ipynb",
        diffexp = "docs/rna-seq-pop-results/notebooks/differential-expression.ipynb",
        qc = "docs/rna-seq-pop-results/notebooks/quality-control.ipynb" if config['QualityControl']['multiqc']['activate'] else [],
        gsea = "docs/rna-seq-pop-results/notebooks/gene-set-enrichment-analysis.ipynb" if config['DifferentialExpression']['GSEA']['activate'] else [],
        gene_families = "docs/rna-seq-pop-results/notebooks/gene-families-heatmap.ipynb" if config['miscellaneous']['GeneFamiliesHeatmap']['activate'] else [],
        pca = "docs/rna-seq-pop-results/notebooks/principal-components-analysis.ipynb" if config['VariantAnalysis']['pca']['activate'] else [],
        genetic_diversity = "docs/rna-seq-pop-results/notebooks/genetic-diversity.ipynb" if config['VariantAnalysis']['summaryStatistics']['activate'] else [],
        selection = "docs/rna-seq-pop-results/notebooks/windowed-selection.ipynb" if config['VariantAnalysis']['selection']['activate'] else [],
        voi = "docs/rna-seq-pop-results/notebooks/variants-of-interest.ipynb" if config['miscellaneous']['VariantsOfInterest']['activate'] else [],
        karyo = "docs/rna-seq-pop-results/notebooks/karyotype.ipynb" if config['VariantAnalysis']['karyotype']['activate'] else [],
    output:
        out_nb = "results/notebooks/process-notebooks.ipynb",
        process = touch("results/rna-seq-pop-results/.processed_nbs")
    conda:
        f'{workflow.basedir}/envs/pythonGenomics.yaml'
    log:
        "logs/notebooks/remove_input_param_cell_nbs.log"
    params:
        wkdir = wkdir
    shell:
        """
        papermill -k pythonGenomics {input.input_nb} {output.out_nb} -p wkdir {params.wkdir} 2> {log}
        """


def getFASTQs(wildcards, rules=None):
    """
    Get FASTQ files from unit sheet.
    If there are more than one wildcard (aka, sample), only return one fastq file
    If the rule is HISAT2align, then return the fastqs with -1 and -2 flags
    """
    metadata = load_metadata(config["metadata"])
    
    if config['fastq']['paired'] == True:
        fastq_cols = ['fq1', 'fq2']
    else:
        fastq_cols = ['fq1']

    if config["QualityControl"]["fastp-trim"]["activate"] == True:
        if rules in ["KallistoQuant", "HISAT2align", "HISAT2align_input"]:
            for i, col in enumerate(fastq_cols):
                metadata = metadata.assign(**{col: f"resources/reads/trimmed/" + metadata["sampleID"] + f"_{i+1}.fastq.gz"})     
            metadata = metadata.set_index("sampleID")
            
            u = metadata.loc[wildcards.sample, fastq_cols].dropna()
            if rules == "HISAT2align":
                return [f"-1 {u.fq1} -2 {u.fq2}"] if config['fastq']['paired'] == True else f"-U {u.fq1}"
            else:
                return [u.fq1, u.fq2] if config['fastq']['paired'] == True else [u.fq1]

    if config["fastq"]["auto"]:
        for i, col in enumerate(fastq_cols):
            metadata = metadata.assign(**{col: f"resources/reads/" + metadata["sampleID"] + f"_{i+1}.fastq.gz"})     
        metadata = metadata.set_index("sampleID")
    else:
        assert (
            "fq1" in metadata.columns
        ), f"The fq1 column in the metadata does not seem to exist. Please create one, or use the 'auto' option and name the fastq files as specified in the config/README.md"
        if config['fastq']['paired']:
            assert (
                "fq2" in metadata.columns
            ), f"The fq2 column in the metadata does not seem to exist. Please create one, or use the 'auto' option and name the fastq files as specified in the config/README.md"
    
        metadata = metadata.set_index("sampleID")

    u = metadata.loc[wildcards.sample, fastq_cols].dropna()
    if rules == "HISAT2align":
        return [f"-1 {u.fq1} -2 {u.fq2}"] if config['fastq']['paired'] == True else f"-U {u.fq1}"
    else:
        return [u.fq1, u.fq2] if config['fastq']['paired'] == True else [u.fq1]

def get_venn_list():
    import itertools

    comp_list = list()
    for all_de_comps in itertools.combinations(config["contrasts"], 2):
        all_de_comps_string = ".".join(all_de_comps)
        comp_list.append(all_de_comps_string)
    for all_de_comps in itertools.combinations(config["contrasts"], 3):
        all_de_comps_string = ".".join(all_de_comps)
        comp_list.append(all_de_comps_string)
    return comp_list


def GetDesiredOutputs(wildcards):

    """
    Function that returns a list of the desired outputs for the rule all, depending on the config.yaml
    configuration file. As of V0.4.0 Does not list every single output, but should mean all rules and desired outputs are created.
    """

    wanted_input = []

    wanted_input.extend(["results/.input.check"])

    # QC & Coverage
    if config["QualityControl"]["multiqc"]['activate']:
        wanted_input.extend(
            expand(
                [
                    "results/qc/multiQC.html",
                ],
                sample=samples,
            )
        )

    if config["QualityControl"]['fastp-trim']["activate"]:
        if config['fastq']['paired'] == True:
            wanted_input.extend(
                expand(
                    [
                        "results/qc/{sample}.html",
                        "results/qc/{sample}.json"
                    ],
                    sample=samples,
                )
            )


    if config["DifferentialExpression"]['gene-level']["activate"]:
        wanted_input.extend(
            expand(
                [
                    "results/genediff/{comp}.csv",
                    "results/genediff/{dataset}_diffexp.xlsx",
                    "results/counts/PCA.pdf",
                    "results/counts/countStatistics.tsv",
                ],
                comp=config["contrasts"],
                dataset=config["dataset"],
            )
        )

    if config["DifferentialExpression"]['isoform-level']["activate"]:
        wanted_input.extend(
            expand(
                [
                    "results/isoformdiff/{comp}.csv",
                    "results/isoformdiff/{dataset}_isoformdiffexp.xlsx",
                ],
                comp=config["contrasts"],
                dataset=config["dataset"],
            )
        )

    if config["DifferentialExpression"]["progressiveGenes"]["activate"]:
        wanted_input.extend(
            expand(
                "results/genediff/{name}.{direction}.progressive.tsv",
                name=config["DifferentialExpression"]["progressiveGenes"]["groups"],
                direction=["up", "down"],
            )
        )

    if config["VariantAnalysis"]["activate"]:
        wanted_input.extend(
            expand(
                [
                    "results/qc/vcfs/{contig}.txt",
                    "results/variantAnalysis/pca/PCA-{contig}-{dataset}.svg",
                    "results/variantAnalysis/SNPstats/snpsPerGenomicFeature.tsv",
                    "results/variantAnalysis/SNPstats/nSNPsPerGene.tsv",
                    "results/variantAnalysis/diversity/{dataset}_SNPdensity_{contig}.svg",
                    "results/variantAnalysis/diversity/SequenceDiversity.tsv",
                    "results/qc/alignments/{sample}.flagstat",
                ],
                contig=config["contigs"],
                dataset=config["dataset"],
                sample=samples,
            )
        )

        if config['QualityControl']['mosdepth']['activate']:
            wanted_input.extend(
                expand(
                    [
                        "results/qc/coverage/{sample}.mosdepth.summary.txt",
                    ],
                sample=samples,
                )
            )

        if config["VariantAnalysis"]["selection"]["activate"]:
            wanted_input.extend(
                expand(
                    [
                        "results/variantAnalysis/selection/FstPerGene.tsv",
                        "results/variantAnalysis/diversity/SequenceDivPerGene.tsv",
                        "results/variantAnalysis/diversity/DxyPerGene.tsv",
                        "results/variantAnalysis/selection/TajimasDPerGene.tsv",
                        "results/variantAnalysis/selection/fst/{wsize}/{comp}.Fst.{contig}.svg",
                    ],
                    contig=config["contigs"],
                    comp=config["contrasts"],
                    wsize=['1000snp_window', '2000snp_window', '5000snp_window'],
                )
            )

    if config["VariantAnalysis"]["ancestry"]["activate"]:
        wanted_input.extend(
            expand(
                [
                    "results/variantAnalysis/ancestry/AIMs_summary.tsv",
                    "results/variantAnalysis/ancestry/AIM_fraction_whole_genome.svg",
                    "results/variantAnalysis/ancestry/n_AIMS_per_chrom.tsv",
                    "results/variantAnalysis/ancestry/AIM_fraction_{contig}.tsv",
                ],
                contig=config["contigs"],
            )
        )

    if config["miscellaneous"]["VariantsOfInterest"]["activate"]:
        wanted_input.extend(
            [
                "results/variantAnalysis/variantsOfInterest/alleleBalance.xlsx",
                "results/variantAnalysis/variantsOfInterest/VOI.heatmapPerSample.svg",
                "results/variantAnalysis/variantsOfInterest/VOI.heatmapPerTreatment.svg",
            ]
        )

    if config["DifferentialExpression"]["GSEA"]["activate"]:
        wanted_input.extend(
            expand(
                [
                    "results/gsea/genediff/{comp}.de.tsv",
                ],
                comp=config["contrasts"],
            )
        )

    if (
        config["VariantAnalysis"]['selection']["activate"]
        and config["DifferentialExpression"]["GSEA"]["activate"]
    ):
        wanted_input.extend(
            expand(
                ["results/gsea/fst/{comp}.fst.tsv"],
                comp=config["contrasts"],
            )
        )

    if config["DifferentialExpression"]["venn"]["activate"]:
        comp_list = get_venn_list()

        wanted_input.extend(
            expand(
                "results/genediff/venn/{comp_string}-{dir_}-Venn.png",
                comp_string=comp_list,
                dir_=["up", "down"],
            )
        )

    if config["VariantAnalysis"]["karyotype"]["activate"]:
        wanted_input.extend(
            expand(
                [
                    "results/karyotype/{karyo}.{dataset}.karyo.txt",
                    "results/karyotype/karyoFreqs.svg",
                ],
                karyo=config["VariantAnalysis"]["karyotype"]["inversions"],
                dataset=config["dataset"],
            )
        )

    if config["miscellaneous"]["sweeps"]["activate"]:
        if config["DifferentialExpression"]['gene-level']["activate"]:
            wanted_input.extend(
                expand(
                    [
                        "results/genediff/ag1000gSweeps/{comp}_swept.tsv",
                    ],
                    comp=config["contrasts"],
                )
            )

    if config["miscellaneous"]["GeneFamiliesHeatmap"]["activate"]:
        wanted_input.extend(
            expand(
                [
                    "results/genediff/GeneFamiliesHeatmap.pdf",
                ],
            )
        )

    if config['jupyter-book']['activate']:
        wanted_input.extend(
            [
                "results/rna-seq-pop-results/_build/html/index.html"
            ]
        )

    return wanted_input


def welcome(version):

    print("---------------------------- RNA-Seq-Pop ----------------------------")
    print(f"Running RNA-Seq-Pop snakemake workflow in {workflow.basedir}\n")
    print(f"Author:   Sanjay Curtis Nagi")
    print(f"Workflow Version: {version}")
    print("Execution time: ", datetime.datetime.now())
    print(f"Dataset: {config['dataset']}", "\n")
