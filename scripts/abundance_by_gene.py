'''
abundance_by_gene.py
Run by sample
- Collapse the Kallisto transcript quantification to gene quantification.
- Calculate totals of intronic and exonic transcripts as summed tpm and summed
estimated counts.
- Recalculate the TPM based on the by gene exon and intron lengths.
- For the filtered version, remove the reads assigned to specified features,
such as spike-in genes and rRNAs, before TPM calculation.
'''

import pandas as pd
import numpy as np
import sys

df = pd.read_csv(snakemake.input['abundance'], sep ='\t')
txt_df = pd.read_csv(snakemake.params['txt_2_gene_file'], sep='\t', header=None, names=['transcript', 'gene', 'symbol', 'blank', 'chrom', 'start', 'end', 'strand'])
len_df = pd.read_csv(snakemake.params['len_file']).rename(columns={'intron_bp':'intron_length', 'exon_bp':'exon_length'})
intron_marker = snakemake.params['intron_marker']
remove_spike_inspect = snakemake.params['remove_spike_inspect']
remove_rrna_inspect = snakemake.params['remove_rrna_inspect']

#Remove spike-in values because they might mess up the inter-library scaling
if remove_spike_inspect == True:
    non_spike = df['target_id'].apply(lambda x: not (x.startswith('SIRV') or x.startswith('ERCC')))
    df = df[non_spike].copy()

#Produce a table of intronic and exonic tpm from the Kallisto output
df['intron'] = df['target_id'].apply(lambda x: intron_marker in x)

# Add the intron and exon lengths into the txt_df
feat_df = pd.merge(txt_df, len_df, left_on='gene', right_on='gene')
df2 = pd.merge(df, feat_df[['transcript', 'gene']], left_on='target_id', right_on='transcript')

# Remove rRNA values because they might mess up the inter-library scaling
if remove_rrna_inspect == True:
    rrna_gene_file = snakemake.params['rrna_gene_file']
    rrna_ids = set(pd.read_csv(rrna_gene_file, header=None)[0].values)
    df2 = df2[~df2.gene.isin(rrna_ids)].copy()

intronic_df = df2[df2['intron']].groupby('gene').agg(intronic_tpm = pd.NamedAgg(column = 'tpm', aggfunc = 'sum'),
                      intronic_est_counts = pd.NamedAgg(column = 'est_counts', aggfunc = 'sum'))
exonic_df = df2[~df2['intron']].groupby('gene').agg(exonic_tpm = pd.NamedAgg(column = 'tpm', aggfunc = 'sum'),
                      exonic_est_counts = pd.NamedAgg(column = 'est_counts', aggfunc = 'sum'))

df3 = pd.concat([intronic_df, exonic_df], axis = 1, sort = False)
df3.fillna(value = 0, inplace = True)
df3 = pd.merge(df3.reset_index(), feat_df[['gene', 'symbol', 'intron_length', 'exon_length']].drop_duplicates(subset=['gene']), left_on='gene', right_on='gene', how='left')
df3['RPL_exon'] = df3['exonic_est_counts']/df3['exon_length']
df3['RPL_intron'] = df3['intronic_est_counts']/df3['intron_length']
#fillna to replace NAs caused by non-existent intron/exon division by 0.
df3[['RPL_exon', 'RPL_intron']] = df3[['RPL_exon', 'RPL_intron']].fillna(value = 0)
RPL_sum = df3['RPL_exon'].sum() + df3['RPL_intron'].sum()
df3['exonic_tpm_recalc'] = 1e6*df3['RPL_exon']/RPL_sum
df3['intronic_tpm_recalc'] = 1e6*df3['RPL_intron']/RPL_sum

#Also combine the intronic and exonic counts
#call it summed_tpm and summed_est_counts to indicate it's summed over all transcripts and intronic, exonic
df3['summed_tpm'] = df3['intronic_tpm'] + df3['exonic_tpm']
df3['summed_tpm_recalc'] = df3['intronic_tpm_recalc'] + df3['exonic_tpm_recalc']
df3['summed_est_counts'] = df3['intronic_est_counts'] + df3['exonic_est_counts']

cols = df3.columns.tolist()
cols.insert(0, cols.pop(cols.index('symbol')))
cols.insert(0, cols.pop(cols.index('gene')))
cols = [i for i in cols if i not in ['RPL_exon', 'RPL_intron']]
df3[cols].to_csv(snakemake.output['gene_table'])