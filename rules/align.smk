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
            return expand('{}/{sample}-{unit}.{group}.fastq.gz'.format(trimmed_dir),
                          group=[1, 2], **wildcards)
        # single end sample
        return '{}/{sample}-{unit}.fastq.gz'.format(trimmed_dir, **wildcards)


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
        star_index = 'indices/star_index_{index_name}'.format(index_name = config['index_name'])
    output:
        'star/{sample}-{unit}/Aligned.out.bam',
        'star/{sample}-{unit}/ReadsPerGene.out.tab'
    log:
        'logs/star/{sample}-{unit}.log'
    params:
        # optional parameters
        extra='--quantMode GeneCounts --sjdbGTFfile {} {}'.format(
              'indices/combo_files/{}.gtf'.format(config['index_name']), config['params']['star'])
    threads: 24
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/run_star.py'

rule quantify_kallisto:
    input:
        fastq = get_fq,
        kallisto_index = 'indices/kallisto_index/{}.idx'.format(config['index_name'])
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

rule summarize_kallisto:
    input:
        abundance = 'kallisto/{sample}-{unit}/abundance.tsv',
        txt_2_gene_file = 'indices/combo_files/{}_txt2gene.txt'.format(config['index_name'])
    output:
        gene_table = 'kallisto/{sample}-{unit}/abundance_by_gene.csv'
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/abundance_by_gene.py'

rule collate_kallisto:
    input:
        infiles = expand('kallisto/{unit.sample}-{unit.unit}/abundance_by_gene.csv', unit=units.itertuples())
    output:
        gene_table = report('results/gene_quantification/summary_abundance_by_gene.csv', '../report/gene_quantification.rst', category = 'Gene Quantification')
    params:
        units_file = config['units']
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/collate_kallisto.py'
