'''
collate_kallisto.py
- Collect the abundance_by_gene values for all experiments and output in a .csv
'''

import pandas as pd

abundance_files = snakemake.input['infiles']
sample_info = snakemake.params['sample_info']

counts = [pd.read_csv(f, index_col = 'gene') for f in abundance_files]

df_list = []
for c, info in zip(counts, sample_info):
    experiment, rep = info
    c['experiment'] = experiment
    c['replicate'] = rep
    df_list.append(c)

#combine all experiments and reps into one file
df = pd.concat(df_list)
df.set_index(['symbol', 'experiment', 'replicate'], append = True, inplace = True)
df = df.reorder_levels(['gene', 'symbol', 'experiment', 'replicate'])
df.to_csv(snakemake.output['gene_table'])
