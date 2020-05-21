# Snakemake workflow: ribopop_rnaseq

This workflow performs analysis of RNA-Seq data for the Ribo-Pop paper.
First sequences are downloaded, and STAR and Kallisto indices are built.
Transcript abundance is estimated with Kallisto and gene-level differential
expression analysis is performed with DESeq2. Reads are also aligned with STAR
and subjected to QC analysis with RSeQC.

The following repository was used as an initial template:
https://github.com/snakemake-workflows/rna-seq-star-deseq2

## Running the workflow

Copy the files from .example/ to the output directory.

Then run the pipeline:

    snakemake --directory <outdir> --use-conda

## Authors

* Mary Kay Thompson (@marykthompson)
