'''
abundance_by_gene.py
- Collapse the Kallisto transcript quantification to gene quantification.
- Calculate totals of primary and mature transcripts as summed tpm and summed
estimated counts.
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

df = pd.read_csv(snakemake.input['abundance'], sep ='\t')

#Produce a table of primary and mature tpm from the Kallisto output
df['intron'] = df['target_id'].apply(lambda x: '-' in x)
df['name'] = df['target_id'].apply(lambda x: x.split('-')[0])
df['gene'] = df['target_id'].map(txt_2_gene)

primary_df = df[df['intron']].groupby('gene').sum()[['tpm', 'est_counts']].rename(columns = {'tpm': 'primary_tpm', 'est_counts': 'primary_est_counts'})
mature_df = df[~df['intron']].groupby('gene').sum()[['tpm', 'est_counts']].rename(columns = {'tpm': 'mature_tpm', 'est_counts': 'mature_est_counts'})
new_df = pd.concat([primary_df, mature_df], axis = 1, sort = False)
new_df.fillna(value = 0, inplace = True)
new_df.index.name = 'gene'
new_df['symbol'] = pd.Series(gene_2_symbol)

#Also combine the primary and mature counts
#call it summed_tpm and summed_est_counts to indicate it's summed over all transcripts and primary, mature
new_df['summed_tpm'] = new_df['primary_tpm'] + new_df['mature_tpm']
new_df['summed_est_counts'] = new_df['primary_est_counts'] + new_df['mature_est_counts']

cols = new_df.columns.tolist()
cols.insert(0, cols.pop(cols.index('symbol')))
new_df[cols].to_csv(snakemake.output['gene_table'])
