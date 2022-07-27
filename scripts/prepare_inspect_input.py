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
nas_exp_des_file = snakemake.output['nas_exp_des_file']
tot_exp_des_file = snakemake.output['tot_exp_des_file']
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

elif analysis_type == 'steadystate':
    df['expname'] = df.apply(lambda x: '%s_%s' % (x['condition'], x['replicate']), axis = 1)
    #https://stackoverflow.com/questions/23482668/sorting-by-a-custom-list-in-pandas
    sorter = condition_mapping
    df.condition = df.condition.astype('category')
    df.condition.cat.set_categories(sorter, inplace=True)
    df.sort_values(['replicate', 'condition'], inplace = True)

nas_df = df[df['RNAtype'] == 'pd'].copy()
tot_df = df[df['RNAtype'] == 'input'].copy()
nas_col_order = nas_df['expname'].unique()
tot_col_order = tot_df['expname'].unique()

if analysis_type == 'timecourse':
    nas_cond = [condition_mapping[int(i.split('_')[0][1:])] for i in nas_col_order]
    tot_cond = [condition_mapping[int(i.split('_')[0][1:])] for i in tot_col_order]
elif analysis_type == 'steadystate':
    nas_cond = [i.split('_')[0] for i in nas_col_order]
    tot_cond = [i.split('_')[0] for i in tot_col_order]

nas_exp_df = pd.DataFrame(nas_cond, columns = ['conditions'])
tot_exp_df = pd.DataFrame(tot_cond, columns = ['conditions'])
nas_exp_df.to_csv(nas_exp_des_file, index = False)
tot_exp_df.to_csv(tot_exp_des_file, index = False)
#write the csv files for INSPEcT
#it writes these as cond1, cond2..., rep2, rep3
#pivot seems to reorder the columns, so need to order them back to the way I specified.
nas_df.reset_index().pivot(index = 'gene', columns = 'expname', values = exonic_col)[nas_col_order].to_csv(nas_exon_file)
nas_df.reset_index().pivot(index = 'gene', columns = 'expname', values = intronic_col)[nas_col_order].to_csv(nas_intron_file)

tot_df.reset_index().pivot(index = 'gene', columns = 'expname', values = exonic_col)[tot_col_order].to_csv(tot_exon_file)
tot_df.reset_index().pivot(index = 'gene', columns = 'expname', values = intronic_col)[tot_col_order].to_csv(tot_intron_file)
