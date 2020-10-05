import pandas as pd

rule DElist:
    output:
        DElist="resources/DE.contrast.list"
    params:
        contrasts = config['contrasts']
    run:
        df = pd.DataFrame(params.contrasts)
        df.columns = ['contrast']
        df.to_csv(output.DElist, sep="\t", index=None)
