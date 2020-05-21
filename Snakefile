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

samples = pd.read_table(config['samples']).set_index('sample', drop = False)
validate(samples, schema = 'schemas/samples.schema.yaml')

units = pd.read_table(config['units'], dtype = str).set_index(['sample', 'unit'], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema = 'schemas/units.schema.yaml')

rnaseq_units = units[units['libtype'] == 'neb_dirII'].copy()

##### target rules #####
rule all:
    input:
        'indices/star_index_{index_name}'.format(index_name = config['index_name']),
        'indices/kallisto_index/{}.idx'.format(config['index_name']),
        expand(['results/diffexp/{contrast}.diffexp.csv',
                'results/diffexp/{contrast}.ma-plot.svg'],
               contrast=config['diffexp']['contrasts']),
        'qc/multiqc_report.html',
        expand('kallisto/{unit.sample}-{unit.unit}/abundance_by_gene.csv', unit = units.itertuples()),
        'results/gene_quantification/summary_abundance_by_gene.csv',
        expand('rrna_coverage/{unit.sample}-{unit.unit}.rrna_bedgraph', unit = rnaseq_units.itertuples())

##### setup report #####
report: 'report/workflow.rst'

##### load rules #####
include: 'rules/prepare_seqs.smk'
include: 'rules/build_indices.smk'
include: 'rules/common.smk'
include: 'rules/trim.smk'
include: 'rules/align.smk'
include: 'rules/diffexp.smk'
include: 'rules/qc.smk'
