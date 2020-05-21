'''
combine_genomes.py
- Combine fasta and gtf files from different genomes, so that can be sent to
indexing software that requires combined files.
'''
import shutil

with open(snakemake.output["combo_fasta"], "wb") as g:
    for i in snakemake.input["fastas"]:
        shutil.copyfileobj(open(i, "rb"), g)
    g.close()

with open(snakemake.output["combo_gtf"], "wb") as g:
    for i in snakemake.input["gtfs"]:
        shutil.copyfileobj(open(i, "rb"), g)
    g.close()
