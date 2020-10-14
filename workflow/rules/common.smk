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


def get_bioc_species_pkg(wildcards):
    """Get the package bioconductor package name for the the species in config.yaml"""
    species_letters = config["resources"]["ref"]["species"][0:2].capitalize()
    return "org.{species}.eg.db".format(species=species_letters)