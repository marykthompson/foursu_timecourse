'''
collate_kallisto.py
- Collect the abundance_by_gene values for all experiments and output in a .csv
'''

import pandas as pd
from collections import defaultdict

units_file = snakemake.params['units_file']
units = pd.read_table(units_file, dtype=str).set_index(['sample', 'unit'], drop=False)
# enforce str in index
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

exps = units.index.levels[0]
reps = units.index.levels[1]

abundance_files = snakemake.input['infiles']
sample_tuples = [(i.sample, i.unit) for i in units.itertuples()]

counts = [pd.read_csv(f, index_col = 'gene') for f in abundance_files]

df_dict = defaultdict(list)
for c, info in zip(counts, sample_tuples):
    experiment, rep = info
    df_dict[experiment].append(c)

exp_dfs = [pd.concat(df_dict[i], keys = reps, names = ['replicate']) for i in exps]
df = pd.concat(exp_dfs, keys = exps, names = ['experiment'])
df.set_index('symbol', append = True, inplace = True)
df = df.reorder_levels(['gene', 'symbol', 'experiment', 'replicate'])
df.to_csv(snakemake.output['gene_table'])
