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

units = pd.read_table(config['units'], dtype = str).set_index(['sample', 'unit'], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema = 'schemas/units.schema.yaml')

##### target rules #####
rule all:
    input:
        'indices/star_index_{index_name}'.format(index_name = config['index_name']),
        'indices/kallisto_index/{}.idx'.format(config['index_name']),
        'qc/multiqc_report.html',
        expand('kallisto/{unit.sample}-{unit.unit}/abundance_by_gene.csv', unit = units.itertuples()),
        'results/gene_quantification/summary_abundance_by_gene.csv',
        'inspect/expDes.csv', 'inspect/nas_exon_tpm.csv', 'inspect/nas_intron_tpm.csv', 'inspect/tot_exon_tpm.csv',
        'inspect/tot_intron_tpm.csv',
        'inspect/synth_rates.csv', 'inspect/deg_rates.csv', 'inspect/proc_rates.csv',
        'inspect/tot_levels.csv', 'inspect/premrna_levels.csv', 'inspect/inspect_data2.rds'

##### setup report #####
report: 'report/workflow.rst'

##### load rules #####
include: 'rules/prepare_seqs.smk'
include: 'rules/build_indices.smk'
include: 'rules/common.smk'
include: 'rules/trim.smk'
include: 'rules/align.smk'
include: 'rules/qc.smk'
include: 'rules/rate_modeling.smk'
