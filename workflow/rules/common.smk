import pandas as pd

rule DElist:
    output:
        "resources/DE.contrast.list"
    params:
        contrasts = config['contrasts']
    run:
        df = pd.DataFrame(snakemake.params[0])
        df.columns = ['contrast']
        df.to_csv(snakemake.output[0], sep="\t", index=None)
