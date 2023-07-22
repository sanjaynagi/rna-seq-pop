rule jupyterbook:
    input:
        pages = "docs/rna-seq-pop-results",
        gsea = "docs/rna-seq-pop-results/notebooks/gene-set-enrichment-analysis.ipynb" if config['DifferentialExpression']['GSEA']['activate'] else [],
        gene_families = "docs/rna-seq-pop-results/notebooks/gene-families-heatmap.ipynb" if config['miscellaneous']['GeneFamiliesHeatmap']['activate'] else [],
        pca = "docs/rna-seq-pop-results/notebooks/principal-components-analysis.ipynb" if config['VariantAnalysis']['pca']['activate'] else [],
        genetic_diversity = "docs/rna-seq-pop-results/notebooks/genetic-diversity.ipynb" if config['VariantAnalysis']['summaryStatistics']['activate'] else [],
        selection = "docs/rna-seq-pop-results/notebooks/windowed-selection.ipynb" if config['VariantAnalysis']['selection']['activate'] else [],
        voi = "docs/rna-seq-pop-results/notebooks/variants-of-interest.ipynb" if config['miscellaneous']['VariantsOfInterest']['activate'] else [],
        karyo = "docs/rna-seq-pop-results/notebooks/karyotype.ipynb" if config['VariantAnalysis']['karyotype']['activate'] else [],
    output:
        directory("results/rna-seq-pop-results/_build/html/"),
        index = "results/rna-seq-pop-results/_build/html/index.html"
    log:
        "logs/jupyterbook/jupyterbook.log"
    conda:
        "../envs/jupyterbook.yaml"
    params: 
        dataset = config['dataset']
    shell:
        """
        jupyter-book build --all {input.pages} --path-output results/rna-seq-pop-results &&
        ln -sf docs/rna-seq-pop-results/_build/html/index.html {params.dataset}-results.html
        """