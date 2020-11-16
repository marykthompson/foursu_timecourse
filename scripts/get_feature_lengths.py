'''
MKT 201015
Parse the txt to gene file and the one of the kallisto quantification files to
produce a dataframes of transcript -> gene and gene -> exon/intron length
'''
import pandas as pd

txt_df = pd.read_csv(snakemake.input['txt_2_gene_file'], sep = '\t',
names = ['transcript_ID', 'gene_ID', 'symbol'], index_col = 'transcript_ID')
#for genes without symbol, assign the gene ID
txt_df['symbol'].fillna(txt_df['gene_ID'], inplace = True)
txt_df.to_pickle(snakemake.output['txt_2_gene_pkl'])

df = pd.read_csv(snakemake.input['kallisto_file'], sep ='\t')
#Produce a table of primary and mature tpm from the Kallisto output
df['intron'] = df['target_id'].apply(lambda x: '-I' in x)
#add in the gene_ID and symbol columns
df = df.merge(txt_df, left_on = 'target_id', right_index = True)
df.rename(columns = {'gene_ID':'gene'}, inplace = True)

#Get length of intronic and exonic regions:
if snakemake.params['feature_length_method'] == 'sum':
    intron_lens = df[df['intron']].groupby('gene')['length'].sum()
    exon_lens = df[~df['intron']].groupby('gene')['length'].sum()
#otherwise use median region lengths
else:
    df['transcript'] = df['target_id'].apply(lambda x: x.split('-I')[0])
    intron_lens = df[df['intron']].groupby(['gene', 'transcript'])['length'].sum().groupby(['gene']).median()
    exon_lens = df[~df['intron']].groupby(['gene', 'transcript'])['length'].sum().groupby(['gene']).median()

exon_lens.name = 'exon_length'
intron_lens.name = 'intron_length'
feat_lens = pd.concat([exon_lens, intron_lens], axis = 1)
#put intron length = 0 for genes with NA
feat_lens.fillna(0, inplace = True)
feat_lens.to_pickle(snakemake.output['feature_len_pkl'])
