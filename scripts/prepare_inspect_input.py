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
excluded_exps = snakemake.params['excluded_exps']

#Load abundance data
df = pd.read_csv(quant_file, index_col = 'gene')

#Remove spike-in values because they might mess up the inter-library scaling
if remove_spike_inspect == True:
    non_spike_idx = df.index.map(lambda x: not (x.startswith('SIRV') or x.startswith('ERCC')))
    df = df.loc[non_spike_idx].copy()

df['exp_type'] = df['experiment'].apply(lambda x: x.split('_')[0])
df['timepoint'] = df['experiment'].apply(lambda x: x.split('_')[1])
df['expname'] = 't' + df['timepoint'] + '_' + df['replicate']
#convert replicate number to integer
df['rep_num'] = df['replicate'].apply(lambda x: x.split('rep')[1])
df['rep_num'] = pd.to_numeric(df['rep_num'], errors = 'coerce')
df['timepoint'] = pd.to_numeric(df['timepoint'])

#drop the experiments that shouldn't be included in the analysis
df = df[~df['experiment'].isin(excluded_exps)].copy()
#drop non-numerical replicates
df.dropna(subset = ['rep_num'], inplace = True)
#sort by replicate and timepoint so that they come out in the correct order
##SHOULDNT THIS BE REPNUM? ISN'T REPLICATE STILL A STRING?
#df.sort_values(by = ['replicate', 'timepoint'], inplace = True)
df.sort_values(by = ['rep_num', 'timepoint'], inplace = True)

#these are the timepoints within one replicate
tpt_array = [tpt_dict[i] for i in df['timepoint'].unique() if i in tpt_dict]
#here also changed from replicate to rep_num
reps = sorted(df['rep_num'].dropna().unique())
num_reps = len(reps)
#write timepoints * replicates for input into R -- this is actually referred to as ExpDf, whereas tpts is without replication.
#Can we use snakemake to give R the tpts and reps directly?
a = tpt_array * num_reps
tdf = pd.DataFrame(a, columns=['timepoints'])
tdf.to_csv(exp_des_file, index=False)
#https://stackoverflow.com/questions/29310792/how-to-save-a-list-as-a-csv-file-with-python-with-new-lines

nas_df = df[df['exp_type'] == 'pd'].copy()
tot_df = df[df['exp_type'] == 'input'].copy()
col_order = df['expname'].unique()

#write the csv files for INSPEcT
#it wrote these as timepoint, timept, rep, rep
#pivot seems to reorder the columns, so need to order them back to the way I specified.
nas_df.reset_index().pivot(index = 'gene', columns = 'expname', values = 'mature_tpm')[col_order].to_csv(nas_exon_file)
nas_df.reset_index().pivot(index = 'gene', columns = 'expname', values = 'primary_tpm')[col_order].to_csv(nas_intron_file)

tot_df.reset_index().pivot(index = 'gene', columns = 'expname', values = 'mature_tpm')[col_order].to_csv(tot_exon_file)
tot_df.reset_index().pivot(index = 'gene', columns = 'expname', values = 'primary_tpm')[col_order].to_csv(tot_intron_file)
