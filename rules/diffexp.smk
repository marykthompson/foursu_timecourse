def get_strandness(units):
    if 'strandedness' in units.columns:
        return units['strandedness'].tolist()
    else:
        strand_list=['none']
        return strand_list*units.shape[0]

def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6

def get_contrast(wildcards):
    return config['diffexp']['contrasts'][wildcards.contrast]

def get_exp_files(wildcards, dir_path = '', file_type = 'abundance.h5'):
    '''
    Get the quantification files for any specific contrast to run.
    Append the path given to the front of the returned files.
    '''
    contrast_samples = get_contrast(wildcards)
    units_w_cond = pd.merge(units, samples[['condition']], left_index = True, right_index = True)
    sample_df = units_w_cond[units_w_cond['condition'].isin(contrast_samples)].copy()
    sample_df['exp_files'] = sample_df.apply(lambda x: os.path.join(dir_path,'-'.join([x['sample'], x['unit']]),file_type), axis = 1)
    return sample_df['exp_files'].tolist()

rule kallisto_deseq2:
    input:
        infiles = lambda wildcards: get_exp_files(wildcards, dir_path = 'kallisto'),
        txt_2_gene_file = 'indices/combo_files/{}_txt2gene.txt'.format(config['index_name'])
    output:
        table = report('results/diffexp/{contrast}.diffexp.csv', '../report/diffexp.rst', category = 'Differential Expression'),
        ma_plot = report('results/diffexp/{contrast}.ma-plot.svg', '../report/ma.rst', category = 'Differential Expression'),
    params:
        contrast = get_contrast,
        units_file = config['units'],
        samples_file = config['samples'],
        #parameterize the sizefactor estimation so that it can be run small or large dataset
        sf_method = config['params']['deseq2']['sf_method']
    conda:
        '../envs/deseq2.yaml'
    log:
        'logs/kallisto_deseq2/{contrast}.diffexp.log'
    threads: 1
    script:
        '../scripts/kallisto_to_deseq2.R'
