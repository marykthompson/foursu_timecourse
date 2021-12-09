def get_fq(wildcards):
    if config['skip_trimming']:
        # no trimming, use raw reads
        return units.loc[(wildcards.sample, wildcards.rep), ['fq1', 'fq2']].dropna()
    else:
        # yes trimming, use trimmed data
        libtype = units.loc[(wildcards.sample, wildcards.rep), 'libtype']
        trimmed_dir = config['params']['trimming'][libtype]['trimmed_dir']
        if not is_single_end(**wildcards):
            return expand('%s/{sample}-{rep}.{group}.fastq.gz' % trimmed_dir,
                           group=[1, 2], **wildcards)
        # single end sample
        return '{}/{sample}-{rep}.fastq.gz'.format(trimmed_dir, **wildcards)

#in order to process the quantseq and rna-seq files separately,
#need to create a list of bb_trimmed/*.fq for the quantseq and cutadapt_trimmed/*.fq
rule align:
    input:
        sample = get_fq,
    output:
        'star/{sample}-{rep}/Aligned.out.bam',
        'star/{sample}-{rep}/ReadsPerGene.out.tab'
    log:
        'logs/star/{sample}-{rep}.log'
    params:
        star_index = config['star_index'],
        gtf = config['gtf_file'],
        options = lambda wildcards: get_program_params(wildcards,
        program = 'star', key = 'options')
    threads: 24
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/run_star.py'

rule quantify_kallisto:
    input:
        fastq = get_fq,
    output:
        'kallisto/{sample}-{rep}/abundance.h5',
        'kallisto/{sample}-{rep}/abundance.tsv',
        'kallisto/{sample}-{rep}/run_info.json'
    log:
        'logs/kallisto/{sample}-{rep}.log'
    params:
        kallisto_index = config['kallisto_index'],
        outdir = 'kallisto/{sample}-{rep}',
        extra = lambda wildcards: get_program_params(wildcards, program = 'kallisto', key = 'options')
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
    input:
        kallisto_file = 'kallisto/{unit.sample}-{unit.replicate}/abundance.tsv'.format(unit = next(units.itertuples())),
        txt_2_gene_file = config['txt_2_gene_file']
    params:
        feature_length_method = config['feature_length_method'],
        feature_length_column = config['feature_length_column']
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
    params:
        remove_spike_inspect = False,
        remove_rrna_inspect = False
    output:
        gene_table = 'kallisto/{sample}-{rep}/abundance_by_gene.csv'
    conda:
        '../envs/main.yaml'
    threads: 8
    script:
        '../scripts/abundance_by_gene.py'

#summarize the kallisto quant after removing specified genes
rule summarize_kallisto_filtered:
    input:
        abundance = 'kallisto/{sample}-{rep}/abundance.tsv',
        txt_2_gene_pkl = 'results/features/txt_2_gene.pkl',
        feature_len_pkl = 'results/features/feature_lens.pkl'
    params:
        remove_spike_inspect = config['remove_spike_inspect'],
        remove_rrna_inspect = config['remove_rrna_inspect'],
        rrna_gene_file = config['rrna_gene_file']
    output:
        gene_table = 'kallisto_filtered/{sample}-{rep}/abundance_by_gene.csv'
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

rule collate_kallisto_filtered:
    input:
        infiles = expand('kallisto_filtered/{unit.sample}-{unit.replicate}/abundance_by_gene.csv', unit=units.itertuples()),
    output:
        gene_table = report('results/gene_quantification/summary_abundance_by_gene_filtered.csv', '../report/gene_quantification_filtered.rst', category = 'Gene Quantification')
    params:
        sample_info = [(i.sample, i.replicate, i.condition, i.RNAtype) for i in units.itertuples()]
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/collate_kallisto.py'
