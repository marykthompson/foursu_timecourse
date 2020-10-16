'''
abundance_by_gene.py
Run by sample
- Collapse the Kallisto transcript quantification to gene quantification.
- Calculate totals of primary and mature transcripts as summed tpm and summed
estimated counts.
- Recalculate the TPM based on the by gene exon and intron lengths.
'''

import pandas as pd
import numpy as np
import sys

df = pd.read_csv(snakemake.input['abundance'], sep ='\t')
txt_df = pd.read_pickle(snakemake.input['txt_2_gene_pkl'])
feat_df = pd.read_pickle(snakemake.input['feature_len_pkl'])

#Produce a table of primary and mature tpm from the Kallisto output
df['intron'] = df['target_id'].apply(lambda x: '_' in x)
df = df.merge(txt_df, left_on = 'target_id', right_index = True)
df.rename(columns = {'gene_ID':'gene'}, inplace = True)

#Get sums of tpm and est_counts by gene
#https://pandas.pydata.org/pandas-docs/stable/user_guide/groupby.html#aggregation
primary_df = df[df['intron']].groupby('gene').agg(primary_tpm = pd.NamedAgg(column = 'tpm', aggfunc = 'sum'),
                      primary_est_counts = pd.NamedAgg(column = 'est_counts', aggfunc = 'sum'))
mature_df = df[~df['intron']].groupby('gene').agg(mature_tpm = pd.NamedAgg(column = 'tpm', aggfunc = 'sum'),
                      mature_est_counts = pd.NamedAgg(column = 'est_counts', aggfunc = 'sum'),
                      symbol = pd.NamedAgg(column = 'symbol', aggfunc = 'first'))

new_df = pd.concat([primary_df, mature_df], axis = 1, sort = False)
new_df.fillna(value = 0, inplace = True)

#recalculate the by gene tpm using the given intron and exon lengths:
new_df = new_df.merge(feat_df, left_index = True, right_index = True)
new_df['RPL_exon'] = new_df['mature_est_counts']/new_df['exon_length']
new_df['RPL_intron'] = new_df['primary_est_counts']/new_df['intron_length']
#fillna to replace NAs caused by non-existent intron/exon zero division by 0.
new_df[['RPL_exon', 'RPL_intron']] = new_df[['RPL_exon', 'RPL_intron']].fillna(value = 0)
RPL_sum = new_df['RPL_exon'].sum() + new_df['RPL_intron'].sum()
new_df['TPM_exon'] = 1e6*new_df['RPL_exon']/RPL_sum
new_df['TPM_intron'] = 1e6*new_df['RPL_intron']/RPL_sum
new_df.index.name = 'gene'

#Also combine the primary and mature counts
#call it summed_tpm and summed_est_counts to indicate it's summed over all transcripts and primary, mature
new_df['summed_tpm'] = new_df['primary_tpm'] + new_df['mature_tpm']
new_df['summed_tpm_recalc'] = new_df['TPM_exon'] + new_df['TPM_intron']
new_df['summed_est_counts'] = new_df['primary_est_counts'] + new_df['mature_est_counts']

cols = new_df.columns.tolist()
cols.insert(0, cols.pop(cols.index('symbol')))
cols = [i for i in cols if i not in ['RPL_exon', 'RPL_intron']]
#why isn't it writing gene in the header?
new_df[cols].to_csv(snakemake.output['gene_table'])
