## gene set enrichment analysis

### This rule and fgsea.R script has been adapted from Johannes Koester and David Laehnemanns Kallisto-Sleuth
### RNA-Seq snakemake workflow https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth

rule fgsea:
    input:
        samples="results/sleuth/samples.tsv",
        diffexp="results/tables/diffexp/{model}.genes-mostsigtrans.diffexp.tsv",
        species_anno=get_bioc_pkg_path,
        gene_sets=config["enrichment"]["fgsea"]["gene_sets_file"]
    output:
        enrichment=report(
            "results/enrichment/fgsea/{model}.all-gene-sets.tsv",
            caption="../report/fgsea-table-all.rst",
            category="Gene set enrichment analysis"
        ),
        rank_ties=report(
            "results/enrichment/fgsea/{model}.rank-ties.tsv",
            caption="../report/fgsea-rank-ties.rst",
            category="Gene set enrichment analysis"
        ),
        significant=report(
            "results/enrichment/fgsea/{model}.sig-gene-sets.tsv",
            caption="../report/fgsea-table-significant.rst",
            category="Gene set enrichment analysis"
        ),
        plot=report(
            "results/enrichment/fgsea/{model}.table-plot.pdf",
            caption="../report/fgsea-table-plot.rst",
            category="Gene set enrichment analysis"
        )
    params:
        bioc_pkg=get_bioc_species_pkg,
        model=get_model,
        gene_set_fdr=config["enrichment"]["fgsea"]["fdr_gene_set"],
        nperm=config["enrichment"]["fgsea"]["nperm"],
        covariate=lambda w: config["diffexp"]["models"][w.model]["primary_variable"]
    conda:
        "../envs/fgsea.yaml"
    log:
        "logs/tables/fgsea/{model}.gene-set-enrichment.log"
    threads: 8
    script:
        "../scripts/fgsea.R"