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
txt_df = pd.read_pickle(snakemake.input['txt_2_gene_pkl'])
feat_df = pd.read_pickle(snakemake.input['feature_len_pkl'])
remove_spike_inspect = snakemake.params['remove_spike_inspect']
remove_rrna_inspect = snakemake.params['remove_rrna_inspect']

#Remove spike-in values because they might mess up the inter-library scaling
if remove_spike_inspect == True:
    non_spike = df['target_id'].apply(lambda x: not (x.startswith('SIRV') or x.startswith('ERCC')))
    df = df[non_spike].copy()

#Produce a table of intronic and exonic tpm from the Kallisto output
df['intron'] = df['target_id'].apply(lambda x: '_I' in x)
df = df.merge(txt_df, left_on = 'target_id', right_index = True)
df.rename(columns = {'gene_ID':'gene'}, inplace = True)

#Remove rRNA values because they might mess up the inter-library scaling
if remove_rrna_inspect == True:
    rrna_gene_file = snakemake.params['rrna_gene_file']
    rrna_ids = set(pd.read_csv(rrna_gene_file, header=None)[0].values)
    df = df[~df.gene.isin(rrna_ids)].copy()

#Get sums of tpm and est_counts by gene
#https://pandas.pydata.org/pandas-docs/stable/user_guide/groupby.html#aggregation
intronic_df = df[df['intron']].groupby('gene').agg(intronic_tpm = pd.NamedAgg(column = 'tpm', aggfunc = 'sum'),
                      intronic_est_counts = pd.NamedAgg(column = 'est_counts', aggfunc = 'sum'))
exonic_df = df[~df['intron']].groupby('gene').agg(exonic_tpm = pd.NamedAgg(column = 'tpm', aggfunc = 'sum'),
                      exonic_est_counts = pd.NamedAgg(column = 'est_counts', aggfunc = 'sum'),
                      symbol = pd.NamedAgg(column = 'symbol', aggfunc = 'first'))

new_df = pd.concat([intronic_df, exonic_df], axis = 1, sort = False)
new_df.fillna(value = 0, inplace = True)

#recalculate the by gene tpm using the given intron and exon lengths:
new_df = new_df.merge(feat_df, left_index = True, right_index = True)
new_df['RPL_exon'] = new_df['exonic_est_counts']/new_df['exon_length']
new_df['RPL_intron'] = new_df['intronic_est_counts']/new_df['intron_length']
#fillna to replace NAs caused by non-existent intron/exon division by 0.
new_df[['RPL_exon', 'RPL_intron']] = new_df[['RPL_exon', 'RPL_intron']].fillna(value = 0)
RPL_sum = new_df['RPL_exon'].sum() + new_df['RPL_intron'].sum()
new_df['exonic_tpm_recalc'] = 1e6*new_df['RPL_exon']/RPL_sum
new_df['intronic_tpm_recalc'] = 1e6*new_df['RPL_intron']/RPL_sum
new_df.index.name = 'gene'

#Also combine the intronic and exonic counts
#call it summed_tpm and summed_est_counts to indicate it's summed over all transcripts and intronic, exonic
new_df['summed_tpm'] = new_df['intronic_tpm'] + new_df['exonic_tpm']
new_df['summed_tpm_recalc'] = new_df['intronic_tpm_recalc'] + new_df['exonic_tpm_recalc']
new_df['summed_est_counts'] = new_df['intronic_est_counts'] + new_df['exonic_est_counts']

cols = new_df.columns.tolist()
cols.insert(0, cols.pop(cols.index('symbol')))
cols = [i for i in cols if i not in ['RPL_exon', 'RPL_intron']]
new_df[cols].to_csv(snakemake.output['gene_table'])
