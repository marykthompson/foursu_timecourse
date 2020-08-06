'''
collate_htseq.py
- Collect the abundance_by_gene values for all experiments and output in a .csv
'''

import pandas as pd

txt_2_gene_file = snakemake.input['txt_2_gene_file']
with open(txt_2_gene_file, 'r') as f:
    txt_2_gene = dict(line.strip().split('\t')[0:2] for line in f)

#gene may not have a symbol, for example if a synthetic spike-in RNA
with open(txt_2_gene_file, 'r') as f:
    gene_2_symbol = {}
    for line in f:
        fields = line.strip().split('\t')
        if len(fields) == 3:
            gene_2_symbol[fields[1]] = fields[2]
        else:
            gene_2_symbol[fields[1]] = fields[1]

abundance_files = snakemake.input['infiles']
sample_info = snakemake.params['sample_info']

counts = [pd.read_csv(f, sep = '\t', names = ['gene', 'counts'], index_col = 'gene') for f in abundance_files]

df_list = []
for c, info in zip(counts, sample_info):
    experiment, rep = info
    c['experiment'] = experiment
    c['replicate'] = rep
    df_list.append(c)

#combine all experiments and reps into one file
df = pd.concat(df_list)
df['symbol'] = df.index.get_level_values('gene').map(gene_2_symbol)
df.set_index(['symbol', 'experiment', 'replicate'], append = True, inplace = True)
df = df.reorder_levels(['gene', 'symbol', 'experiment', 'replicate'])
df.to_csv(snakemake.output['gene_table'])
