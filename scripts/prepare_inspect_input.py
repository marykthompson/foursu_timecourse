#prepare_inspect_input
#convert the Kallisto gene quantifications to input tables for the INSPEcT package
import os
import pandas as pd
import numpy as np
import sys

tpt_dict = snakemake.params['tpts']
quant_file = snakemake.input['gene_quant_file']
exp_des_file = snakemake.output['exp_des_file']
nas_exon_file = snakemake.output['nas_exon_file']
nas_intron_file = snakemake.output['nas_intron_file']
tot_exon_file = snakemake.output['tot_exon_file']
tot_intron_file = snakemake.output['tot_intron_file']
remove_spike_inspect = snakemake.params['remove_spike_inspect']

#Load abundance data
df = pd.read_csv(quant_file, index_col = 'gene')

#Remove spike-in values because they might mess up the inter-library scaling
if remove_spike_inspect == True:
    non_spike_idx = df.index.map(lambda x: not (x.startswith('SIRV') or x.startswith('ERCC')))
    df = df.loc[non_spike_idx].copy()

df['exp_type'] = df['experiment'].apply(lambda x: x.split('_')[0])
df['timepoint'] = df['experiment'].apply(lambda x: x.split('_')[1])
df['expname'] = 't' + df['timepoint'] + '_' + df['replicate']

#simply for the test, ideally use integers in the units.tsv file in the pipeline and then don't need:
df['replicate'].replace(to_replace = 'rep1', value = '1', inplace = True)
#convert replicate and timepoint to numeric values and sort.
df['replicate'] = pd.to_numeric(df['replicate'], errors = 'coerce')
df['timepoint'] = pd.to_numeric(df['timepoint'])

#drop the input5 sample because it doesn't have a matched pulldown sample
df = df[df['experiment'] !=  'input_5'].copy()

#drop the replicate 1a one
df.dropna(subset = ['replicate'], inplace = True)

#copy df to make a df2 which would be a like a separate replicate
df2 = df.copy()
#in the real experiment df2 would be part of the same df, so concatenate it back at the end.

#how can we stimulate noise for the second replicate?
#https://stackoverflow.com/questions/14058340/adding-noise-to-a-signal-in-python
#https://stackoverflow.com/questions/46093073/adding-gaussian-noise-to-a-dataset-of-floating-points-and-save-it-python?rq=1

mu, sigma = 0, 1
noise1 = np.random.normal(mu, sigma, len(df2))
noise2 = np.random.normal(mu, sigma, len(df2))

df2['primary_tpm'] = np.clip(df2['primary_tpm'] + noise1, 0, None)
df2['mature_tpm'] = np.clip(df2['mature_tpm'] + noise2, 0, None)

df2['replicate'] = 2.0
cdf = pd.concat([df, df2])
#repeat the expname assignment
cdf['expname'] = cdf.apply(lambda x: 't%s_%s' % (x['timepoint'], int(x['replicate'])), axis = 1)
#sort by replicate and timepoint so that they come out in the correct order
cdf.sort_values(by = ['replicate', 'timepoint'], inplace = True)

#these are the timepoints within one replicate
tpt_array = [tpt_dict[i] for i in cdf['timepoint'].unique() if i in tpt_dict]
reps = sorted(cdf['replicate'].dropna().unique())
num_reps = len(reps)
#write timepoints * replicates for input into R -- this is actually referred to as ExpDf, whereas tpts is without replication.
#Can we use snakemake to give R the tpts and reps directly?
a = tpt_array * num_reps
tdf = pd.DataFrame(a, columns=['timepoints'])
tdf.to_csv(exp_des_file, index=False)
#https://stackoverflow.com/questions/29310792/how-to-save-a-list-as-a-csv-file-with-python-with-new-lines

nas_df = cdf[cdf['exp_type'] == 'pd'].copy()
tot_df = cdf[cdf['exp_type'] == 'input'].copy()
col_order = cdf['expname'].unique()

#write the csv files for INSPEcT
#it wrote these as timepoint, timept, rep, rep
#pivot seems to reorder the columns, so need to order them back to the way I specified.
nas_df.reset_index().pivot(index = 'gene', columns = 'expname', values = 'mature_tpm')[col_order].to_csv(nas_exon_file)
nas_df.reset_index().pivot(index = 'gene', columns = 'expname', values = 'primary_tpm')[col_order].to_csv(nas_intron_file)

tot_df.reset_index().pivot(index = 'gene', columns = 'expname', values = 'mature_tpm')[col_order].to_csv(tot_exon_file)
tot_df.reset_index().pivot(index = 'gene', columns = 'expname', values = 'primary_tpm')[col_order].to_csv(tot_intron_file)
