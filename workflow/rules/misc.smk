
comp_list = get_venn_list()

rule VennDiagrams:
    """
    Find intersection of DE analyses between comparisons and plot
    """
    input:
        DE=expand("results/genediff/{comp}.csv", comp=config["contrasts"]),
    output:
        venn=expand(
            "results/genediff/venn/{comparisons}-{{dir_}}-Venn.png",
            comparisons=comp_list,
        ),
    conda:
        "../envs/venn.yaml"
    log:
        "logs/VennDiagrams_{dir_}.log",
    params:
        comparisons=config["contrasts"],
        padj_threshold=config["DifferentialExpression"]["venn"]["padj_threshold"],
        dataset=config["dataset"],
    script:
        "../scripts/VennDiagrams.py"
