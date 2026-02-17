if config["VariantAnalysis"]["selection"]["population-branch-statistic"]["activate"]:
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



def get_fastqs(wildcards, rules=None):
    """
    Get FASTQ files from unit sheet.
    If there are more than one wildcard (aka, sample), only return one fastq file
    If the rule is hisat2_align, then return the fastqs with -1 and -2 flags
    """
    metadata = load_metadata(config["metadata"])
    
    if config['fastq']['paired'] == True:
        fastq_cols = ['fq1', 'fq2']
    else:
        fastq_cols = ['fq1']

    if config["QualityControl"]["fastp-trim"]["activate"] == True:
        if rules in ["kallisto_quant", "hisat2_align", "hisat2_align_input"]:
            for i, col in enumerate(fastq_cols):
                metadata = metadata.assign(**{col: f"resources/reads/trimmed/" + metadata["sampleID"] + f"_{i+1}.fastq.gz"})     
            metadata = metadata.set_index("sampleID")
            
            u = metadata.loc[wildcards.sample, fastq_cols].dropna()
            if rules == "hisat2_align":
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
    if rules == "hisat2_align":
        return [f"-1 {u.fq1} -2 {u.fq2}"] if config['fastq']['paired'] == True else f"-U {u.fq1}"
    else:
        return [u.fq1, u.fq2] if config['fastq']['paired'] == True else [u.fq1]


def get_bam(wildcards):
    """
    Get BAM files depending on aligner
    """
    if config['pipeline'] == 'parabricks':
        bam = "results/alignments/{sample}.star.bam"
    else:
        bam = "results/alignments/{sample}.hisat2.bam"
    return bam

def rnaseqpop_outputs(wildcards):

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

        if config["VariantAnalysis"]["karyotype"]["activate"]:
            wanted_input.extend(
                expand(
                    [
                        "results/karyotype/karyotypes.tsv",
                        "results/karyotype/karyotype_heatmap.png",
                    ],
                    karyo=config["VariantAnalysis"]["karyotype"]["inversions"],
                    dataset=config["dataset"],
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

    if config['QualityControl']['coverage']['activate']:
        wanted_input.extend(
            expand(
                [
                    "results/qc/coverage/{sample}.mosdepth.summary.txt",
                ],
            sample=samples,
            )
        )


    if config["VariantsOfInterest"]["activate"]:
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

    if config['results-jupyterbook']['activate']:
        wanted_input.extend(
            [
                "results/rna-seq-pop-results/_build/html/index.html"
            ]
        )

    return wanted_input


def welcome(version, using_user_metadata_colours=False):
    import datetime

    try:
        from rich.console import Console
        from rich.panel import Panel
        from rich.table import Table
        use_rich = True
    except ImportError:
        use_rich = False

    now = datetime.datetime.now().replace(microsecond=0)

    fastq_cfg = config.get("fastq", {})
    paired = fastq_cfg.get("paired", True)
    auto_fastq = fastq_cfg.get("auto", False)
    qc_cfg = config.get("QualityControl", {})
    trim_active = qc_cfg.get("fastp-trim", {}).get("activate", False)
    pipeline_mode = config.get("pipeline", "cpu")

    if trim_active:
        input_str = "FASTQs from resources/reads/trimmed/ generated by fastp-trim"
    elif auto_fastq and paired:
        input_str = "Paired-end FASTQs auto-detected in resources/reads/ (<sampleID>_1/_2.fastq.gz)"
    elif auto_fastq and not paired:
        input_str = "Single-end FASTQs auto-detected in resources/reads/ (<sampleID>_1.fastq.gz)"
    elif paired:
        input_str = "Paired-end FASTQ paths provided in metadata columns fq1/fq2"
    else:
        input_str = "Single-end FASTQ paths provided in metadata column fq1"

    aligner_str = "STAR + GATK HaplotypeCaller (Parabricks)" if pipeline_mode == "parabricks" else "HISAT2 + FreeBayes"

    active_modules = []
    if qc_cfg.get("multiqc", {}).get("activate", False):
        active_modules.append("MultiQC")
    if qc_cfg.get("coverage", {}).get("activate", False):
        active_modules.append("Coverage")
    if config.get("DifferentialExpression", {}).get("gene-level", {}).get("activate", False):
        active_modules.append("DE genes")
    if config.get("DifferentialExpression", {}).get("isoform-level", {}).get("activate", False):
        active_modules.append("DE isoforms")
    if config.get("DifferentialExpression", {}).get("GSEA", {}).get("activate", False):
        active_modules.append("GSEA")
    if config.get("VariantAnalysis", {}).get("activate", False):
        active_modules.append("Variant analysis")
    if config.get("VariantsOfInterest", {}).get("activate", False):
        active_modules.append("Variants of interest")
    if config.get("results-jupyterbook", {}).get("activate", False):
        active_modules.append("Results jupyter-book")

    active_modules_str = ", ".join(active_modules) if active_modules else "None"
    contigs = config.get("contigs", [])
    contrasts = config.get("contrasts", [])
    sample_count = len(samples)
    voi_path = config.get("VariantsOfInterest", {}).get("path", "N/A")

    if use_rich:
        console = Console()
        table = Table.grid(padding=(0, 2))
        table.add_row("Working dir:", workflow.basedir)
        table.add_row("Author:", "Sanjay Curtis Nagi")
        table.add_row("Workflow version:", version)
        table.add_row("Execution time:", str(now))
        table.add_row("Dataset:", str(config.get("dataset", "N/A")))
        table.add_row("Metadata file:", str(config.get("metadata", "N/A")))
        table.add_row("Samples:", str(sample_count))
        table.add_row("Pipeline:", aligner_str)
        table.add_row("Input:", input_str)
        table.add_row("Contigs:", ", ".join(map(str, contigs)) if contigs else "N/A")
        table.add_row("Contrasts:", ", ".join(map(str, contrasts)) if contrasts else "N/A")
        table.add_row("Active modules:", active_modules_str)
        if config.get("VariantsOfInterest", {}).get("activate", False):
            table.add_row("VOI table:", str(voi_path))
        if using_user_metadata_colours:
            table.add_row("Metadata colours:", "user-provided schema in config/metadata_colours.json")

        console.print(
            Panel(
                table,
                title="[bold sea_green2]RNA-Seq-Pop[/bold sea_green2]",
                border_style="sea_green2",
            )
        )
    else:
        sep = "-" * 78
        print(f"\n{sep}")
        print(" RNA-Seq-Pop ")
        print(f"{sep}")
        print(f"{'Working dir:':18} {workflow.basedir}")
        print(f"{'Author:':18} Sanjay Curtis Nagi")
        print(f"{'Workflow version:':18} {version}")
        print(f"{'Execution time:':18} {now}")
        print(f"{'Dataset:':18} {config.get('dataset', 'N/A')}")
        print(f"{'Metadata file:':18} {config.get('metadata', 'N/A')}")
        print(f"{'Samples:':18} {sample_count}")
        print(f"{'Pipeline:':18} {aligner_str}")
        print(f"{'Input:':18} {input_str}")
        print(f"{'Contigs:':18} {', '.join(map(str, contigs)) if contigs else 'N/A'}")
        print(f"{'Contrasts:':18} {', '.join(map(str, contrasts)) if contrasts else 'N/A'}")
        print(f"{'Active modules:':18} {active_modules_str}")
        if config.get("VariantsOfInterest", {}).get("activate", False):
            print(f"{'VOI table:':18} {voi_path}")
        if using_user_metadata_colours:
            print(f"{'Metadata colours:':18} user-provided schema in config/metadata_colours.json")
        print(f"{sep}\n")
