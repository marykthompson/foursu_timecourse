def get_fq(wildcards):
    if config['trimming']['skip']:
        # no trimming, use raw reads
        return units.loc[(wildcards.sample, wildcards.unit), ['fq1', 'fq2']].dropna()
    else:
        # yes trimming, use trimmed data
        libtype = units.loc[(wildcards.sample, wildcards.unit), 'libtype']
        trimmed_dir = config['trimming'][libtype]['trimmed_dir']

        if not is_single_end(**wildcards):
            # paired-end sample
            return expand("{trimmed}/{sample}-{unit}.{group}.fastq.gz",
                          group=[1, 2], trimmed = trimmed_dir, **wildcards)

        # single end sample
        return '{trimmed}/{sample}-{unit}.fastq.gz'.format(trimmed = trimmed_dir, **wildcards)


def get_program_params(wildcards, program = ''):
    '''
    Get the params by libtype. Different libtypes will have different
    mapping strands, etc.
    '''
    libtype = units.loc[(wildcards.sample, wildcards.unit), 'libtype']
    extra = config['params'][program][libtype]
    return extra

#in order to process the quantseq and rna-seq files separately,
#need to create a list of bb_trimmed/*.fq for the quantseq and cutadapt_trimmed/*.fq
rule align:
    input:
        sample = get_fq,
        star_index = config['star_index']
    output:
        'star/{sample}-{unit}/Aligned.out.bam',
        'star/{sample}-{unit}/ReadsPerGene.out.tab'
    log:
        'logs/star/{sample}-{unit}.log'
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
        'kallisto/{sample}-{unit}/abundance.h5',
        'kallisto/{sample}-{unit}/abundance.tsv',
        'kallisto/{sample}-{unit}/run_info.json'
    log:
        'logs/kallisto/{sample}-{unit}.log'
    params:
        outdir = 'kallisto/{sample}-{unit}',
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
        kallisto_file = 'kallisto/{unit.sample}-{unit.unit}/abundance.tsv'.format(unit = next(units.itertuples())),
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
        abundance = 'kallisto/{sample}-{unit}/abundance.tsv',
        txt_2_gene_pkl = 'results/features/txt_2_gene.pkl',
        feature_len_pkl = 'results/features/feature_lens.pkl'
    output:
        gene_table = 'kallisto/{sample}-{unit}/abundance_by_gene.csv'
    conda:
        '../envs/main.yaml'
    threads: 8
    script:
        '../scripts/abundance_by_gene.py'

rule collate_kallisto:
    input:
        infiles = expand('kallisto/{unit.sample}-{unit.unit}/abundance_by_gene.csv', unit=units.itertuples()),
    output:
        gene_table = report('results/gene_quantification/summary_abundance_by_gene.csv', '../report/gene_quantification.rst', category = 'Gene Quantification')
    params:
        sample_info = [(i.sample, i.unit) for i in units.itertuples()]
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/collate_kallisto.py'

rule htseq_count:
    input:
        'star/{sample}-{unit}/Aligned.out.bam'
    output:
        'htseq/{sample}-{unit}/htseq_count.txt'
    params:
        combo_gtf = config['gtf_file']
    conda:
        '../envs/main.yaml'
    log:
        'logs/htseq/{sample}-{unit}.log'
    shell:
        'htseq-count {input} {params.combo_gtf} > {output}'

rule collate_htseq:
    input:
        infiles = expand('htseq/{unit.sample}-{unit.unit}/htseq_count.txt', unit=units.itertuples()),
        txt_2_gene_file = config['txt_2_gene_file']
    output:
        gene_table = report('results/gene_quantification/summary_abundance_by_gene_htseq.csv', '../report/gene_quantification_htseq.rst', category = 'Gene Quantification')
    params:
        sample_info = [(i.sample, i.unit) for i in units.itertuples()]
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/collate_htseq.py'
