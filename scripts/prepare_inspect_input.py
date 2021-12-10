'''
Prepare_inspect_input
Convert the Kallisto gene quantifications to input tables for the INSPEcT package.
Use the columns given by intronic_col and exonic_col to choose which metric to output.
'''
import os
import pandas as pd
import numpy as np
import sys

condition_mapping = snakemake.params['condition_mapping']
analysis_type = snakemake.params['analysis_type']
quant_file = snakemake.input['gene_quant_file']
exp_des_file = snakemake.output['exp_des_file']
nas_exon_file = snakemake.output['nas_exon_file']
nas_intron_file = snakemake.output['nas_intron_file']
tot_exon_file = snakemake.output['tot_exon_file']
tot_intron_file = snakemake.output['tot_intron_file']
excluded_exps = snakemake.params['excluded_exps']
intronic_col = snakemake.params['intronic_col']
exonic_col = snakemake.params['exonic_col']

#Load abundance data
df = pd.read_csv(quant_file, index_col = 'gene')

#drop the experiments that shouldn't be included in the analysis
df = df[~df['sample'].isin(excluded_exps)].copy()

#enforce numeric replicate numbers and sort by replicate
df['replicate'] = pd.to_numeric(df['replicate'], errors = 'coerce')
df.dropna(subset = ['replicate'], inplace = True)
num_reps = len(df['replicate'].unique())

if analysis_type == 'timecourse':
    df['expname'] = df.apply(lambda x: 't%s_%s' % (x['condition'], x['replicate']), axis = 1)
    df['condition'] = pd.to_numeric(df['condition'], errors = 'coerce')
    df.dropna(subset = ['condition'], inplace = True)
    #sorting needs to be done at the same time:
    df.sort_values(by = ['replicate', 'condition'], inplace = True)
    #Write the expDes file
    #these are the timepoints within one replicate
    tpt_array = [condition_mapping[i] for i in df['condition'].unique() if i in condition_mapping]
    #write timepoints * replicates for input into R
    #this is actually referred to as ExpDf, whereas tpts is without replication.
    tdf = pd.DataFrame(tpt_array * num_reps, columns=['conditions'])
    tdf.to_csv(exp_des_file, index=False)

elif analysis_type == 'steadystate':
    df['expname'] = df.apply(lambda x: '%s_%s' % (x['condition'], x['replicate']), axis = 1)
    #https://stackoverflow.com/questions/23482668/sorting-by-a-custom-list-in-pandas
    sorter = condition_mapping
    df.condition = df.condition.astype('category')
    df.condition.cat.set_categories(sorter, inplace=True)
    df.sort_values(['replicate', 'condition'], inplace = True)
    #INSPEcT just needs the ordered names of the conditions
    cdf = pd.DataFrame(condition_mapping*num_reps, columns = ['conditions'])
    cdf.to_csv(exp_des_file, index = False)

nas_df = df[df['RNAtype'] == 'pd'].copy()
tot_df = df[df['RNAtype'] == 'input'].copy()
col_order = df['expname'].unique()

#write the csv files for INSPEcT
#it writes these as cond1, cond2..., rep2, rep3
#pivot seems to reorder the columns, so need to order them back to the way I specified.
nas_df.reset_index().pivot(index = 'gene', columns = 'expname', values = exonic_col)[col_order].to_csv(nas_exon_file)
nas_df.reset_index().pivot(index = 'gene', columns = 'expname', values = intronic_col)[col_order].to_csv(nas_intron_file)

tot_df.reset_index().pivot(index = 'gene', columns = 'expname', values = exonic_col)[col_order].to_csv(tot_exon_file)
tot_df.reset_index().pivot(index = 'gene', columns = 'expname', values = intronic_col)[col_order].to_csv(tot_intron_file)
