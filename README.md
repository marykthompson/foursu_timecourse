# Snakemake workflow: foursu_timecourse

This workflow performs analysis of RNA-Seq data from a 4sU treatment time course,
where foursu is applied over a window at different points in a stress response.
First sequences are downloaded, and STAR and Kallisto indices are built.
Transcript abundance is estimated with Kallisto and then collated to produce gene-level
abundance estimates. Abundance estimates are fed into the INSPEcT R package to
estimate synthesis, splicing, and degradation rates. Reads are also aligned with STAR
and subjected to QC analysis with RSeQC.

The following repository was used as an initial template:
https://github.com/snakemake-workflows/rna-seq-star-deseq2

## Running the workflow

Copy the files from .example/ to the output directory.

Then run the pipeline:

    snakemake --directory <outdir> --use-conda

## Authors

* Mary Kay Thompson (@marykthompson)
