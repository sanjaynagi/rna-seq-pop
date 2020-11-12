import pandas as pd

#def check_chroms(gffpath, chroms, str_remove):
#            """
#            Check that the chroms in the gff match ours, at least after removing the chosen string (str_remove)
#            """
#            
#            gff = allel.gff3_to_dataframe(gffpath)
# #           gff['seqid'] = gff['seqid'].str.strip(str_remove)
#                
#            assert np.isin(chroms, gff.seqid).all(), f"All provided chromosome names ({chroms}) are not present in GFF file"
#    

#rule check_inputs:
#    output:
##       touch("config/.input.check")
#    params:
# #      metadata=config['samples'],
#  #     chroms=config['chroms'],
#       gffpath=config['ref']['gff'],
#   ##    str_remove=config['ref']['str_remove'],
#       gene_names=config['ref']['genenames']
#    run:
#        metadata = pd.read_csv(params.metadata, sep="\t")
#        
#        #check to see if fastq files 
#        for sample in metadata.samples:
#            for n in [1,2]:
#                fqpath = f"resources/reads/{sample}_{n}.fastq.gz"
##
#                assert os.path.isfile(fqpath), "all sample names in 'samples.tsv' do not match a .fastq.gz file"#
#
# #       #check that the chromosomes are present in the gff file 
#  #      check_chroms(params.gffpath, params.chroms, params.str_remove)
#
# ##       #check column names of gene_names.tsv
#   #     gene_names = pd.read_csv(params.gene_names, sep="\t")
#    #    colnames = ['Gene_stable_ID', 'Gene_description', 'Gene_name']
#    #    assert (gene_names.columns == colnames).all(), f"Column names of gene_names.tsv should be {colnames}"


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
