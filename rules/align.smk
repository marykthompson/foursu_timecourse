def get_fq(wildcards):
    if config['trimming']['skip']:
        # no trimming, use raw reads
        return units.loc[(wildcards.sample, wildcards.rep), ['fq1', 'fq2']].dropna()
    else:
        # yes trimming, use trimmed data
        libtype = units.loc[(wildcards.sample, wildcards.rep), 'libtype']
        trimmed_dir = config['trimming'][libtype]['trimmed_dir']

        if not is_single_end(**wildcards):
            # paired-end sample
            return expand("{trimmed}/{sample}-{rep}.{group}.fastq.gz",
                          group=[1, 2], trimmed = trimmed_dir, **wildcards)

        # single end sample
        return '{trimmed}/{sample}-{rep}.fastq.gz'.format(trimmed = trimmed_dir, **wildcards)


def get_program_params(wildcards, program = ''):
    '''
    Get the params by libtype. Different libtypes will have different
    mapping strands, etc.
    '''
    libtype = units.loc[(wildcards.sample, wildcards.rep), 'libtype']
    extra = config['params'][program][libtype]
    return extra

#in order to process the quantseq and rna-seq files separately,
#need to create a list of bb_trimmed/*.fq for the quantseq and cutadapt_trimmed/*.fq
rule align:
    input:
        sample = get_fq,
        star_index = config['star_index']
    output:
        'star/{sample}-{rep}/Aligned.out.bam',
        'star/{sample}-{rep}/ReadsPerGene.out.tab'
    log:
        'logs/star/{sample}-{rep}.log'
    params:
        # optional parameters
        extra='--quantMode GeneCounts --sjdbGTFfile {} {}'.format(
              config['gtf_file'], config['params']['star'])
    threads: 24
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/run_star.py'

rule quantify_kallisto:
    input:
        fastq = get_fq,
        kallisto_index = config['kallisto_index']
    output:
        'kallisto/{sample}-{rep}/abundance.h5',
        'kallisto/{sample}-{rep}/abundance.tsv',
        'kallisto/{sample}-{rep}/run_info.json'
    log:
        'logs/kallisto/{sample}-{rep}.log'
    params:
        outdir = 'kallisto/{sample}-{rep}',
        extra = lambda wildcards: get_program_params(wildcards, program = 'kallisto')
    threads: 1
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/run_kallisto_quant.py'

rule get_feature_lengths:
    '''
    Parse the Kallisto output and generate a table with the
    by gene intron and exon lengths.
    We only need to read one of the quantification files because
    they all use the same genome.
    Also parse the txt_2_gene mapping file and store as a pkl.
    '''
    #Do you need to reparse the units file or will it still be available here?
    input:
        kallisto_file = 'kallisto/{unit.sample}-{unit.replicate}/abundance.tsv'.format(unit = next(units.itertuples())),
        txt_2_gene_file = config['txt_2_gene_file']
    params:
        feature_length_method = 'median'
    output:
        txt_2_gene_pkl = 'results/features/txt_2_gene.pkl',
        feature_len_pkl = 'results/features/feature_lens.pkl'
    script:
        '../scripts/get_feature_lengths.py'

rule summarize_kallisto:
    input:
        abundance = 'kallisto/{sample}-{rep}/abundance.tsv',
        txt_2_gene_pkl = 'results/features/txt_2_gene.pkl',
        feature_len_pkl = 'results/features/feature_lens.pkl'
    output:
        gene_table = 'kallisto/{sample}-{rep}/abundance_by_gene.csv'
    conda:
        '../envs/main.yaml'
    threads: 8
    script:
        '../scripts/abundance_by_gene.py'

rule collate_kallisto:
    input:
        infiles = expand('kallisto/{unit.sample}-{unit.replicate}/abundance_by_gene.csv', unit=units.itertuples()),
    output:
        gene_table = report('results/gene_quantification/summary_abundance_by_gene.csv', '../report/gene_quantification.rst', category = 'Gene Quantification')
    params:
        sample_info = [(i.sample, i.replicate, i.condition, i.RNAtype) for i in units.itertuples()]
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/collate_kallisto.py'
