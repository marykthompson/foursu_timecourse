import pandas as pd
import os
from snakemake.utils import validate, min_version

#need to use workflow.basedir to get the relative directory
snake_dir = workflow.basedir

##### set minimum snakemake version #####
min_version('5.1.2')

##### load config and sample sheets #####
configfile: 'config.yaml'
validate(config, schema = 'schemas/config.schema.yaml')

units = pd.read_csv(config['units'], dtype = str).set_index(['sample', 'replicate'], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema = 'schemas/units.schema.yaml')

##### target rules #####
rule all:
    input:
        # 'qc/multiqc_report.html',
        # expand('kallisto/{unit.sample}-{unit.unit}/abundance_by_gene.csv', unit = units.itertuples()),
        # 'results/gene_quantification/summary_abundance_by_gene.csv',
        # 'inspect/nas_exon_tpm.csv', 'inspect/nas_intron_tpm.csv',
         'inspect/synth_rates.csv', 'inspect/deg_rates.csv', 'inspect/proc_rates.csv',
        # list(set(expand('inspect/inspect_data1_{unit.replicate}.rds', unit = units.itertuples()))),
        # 'inspect/inspect_data2.rds',
        # 'results/gene_quantification/summary_abundance_by_gene.csv',
        #list(set(expand('dexseq/{unit.sample}-{unit.replicate}.txt', unit = units.itertuples())))
        # expand('kallisto/{unit.sample}-{unit.replicate}/abundance_by_gene.csv', unit=units.itertuples()),
        # expand('bigwig/{unit.sample}-{unit.replicate}.{strand}.bw', unit = units.itertuples(), strand = ['p', 'm'])

##### setup report #####
report: 'report/workflow.rst'

##### load rules #####
include: 'rules/common.smk'
include: 'rules/trim.smk'
include: 'rules/align.smk'
include: 'rules/qc.smk'
include: 'rules/rate_modeling.smk'
include: 'rules/genome_plots.smk'
include: 'rules/dexseq.smk'
