# Snakemake workflow: foursu_timecourse

This workflow performs analysis of RNA-Seq data from a 4sU treatment time course,
where foursu is applied over a window at different points after a stimulus.
It can also be used to examine kinetics from steady-state experiments.
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

## Constructing the experiment template file (units.csv)

See the example units.csv file. The column definitions are:
- sample*: a name describing the sample.
- replicate*: the replicate, written as a number, e.g. 1
- condition: if a timepoint, this should be a number, 1 to n. For steady-state
comparison it should be a string (e.g. wt, mutant).
- RNAtype: input or pd (pulldown of 4sU-labeled RNA)
- libtype: the name of a library kit, with associated parameters defined in the config file.
- fq1: path to fasta of the mate1 read.
- fq2: path to fasta of the mate2 read, if using paired-end data. For single-end
data, leave the values of this column empty.

*Note that Snakemake will use sample and replicate to create the output directory
for quantitation, so the combination of the two should yield a unique name in the
dataset. A simple way to name the samples is by RNAtype_condition (e.g. input_wt/input_1)
which would result in a folder name for replicate 1 as input_wt-1 or input_1-1.

## Authors

* Mary Kay Thompson (@marykthompson)
