def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.rep), ['fq1', 'fq2']].dropna()

def adapter_string_pe(wildcards, program = '', a = 'adapter1_3p', A = 'adapter2_3p'):
    '''Make adapter string needed for cutadapt-pe'''
    return '-a {} -A {}'.format(get_program_params(wildcards, program = program, key = a),
    get_program_params(wildcards, program = program, key = A))

def adapter_string_se(wildcards, program = '', a = 'adapter1_3p', options = 'options'):
    '''Make adapter string needed for cutadapt-se'''
    return '-a {} {}'.format(get_program_params(wildcards, program = program, key = a),
    get_program_params(wildcards, program = program, key = 'options'))

rule cutadapt_pe:
    input:
        get_fastq
    output:
        fastq1 = 'trimmed/{sample}-{rep}.1.fastq.gz',
        fastq2 = 'trimmed/{sample}-{rep}.2.fastq.gz',
        qc = 'trimmed/{sample}-{rep}.qc.txt'
    params:
        adapters = lambda wildcards: adapter_string_pe(wildcards, program = 'cutadapt'),
        others = lambda wildcards: get_program_params(wildcards, program = 'cutadapt',
        key = 'options')
    log:
        'logs/cutadapt/{sample}-{rep}.log'
    wrapper:
        '0.49.0/bio/cutadapt/pe'

rule cutadapt:
    input:
        get_fastq
    output:
        fastq = 'trimmed/{sample}-{rep}.fastq.gz',
        qc = 'trimmed/{sample}-{rep}.qc.txt'
    params:
        lambda wildcards: adapter_string_se(wildcards, program = 'cutadapt')
    log:
        'logs/cutadapt/{sample}-{rep}.log'
    wrapper:
        '0.49.0/bio/cutadapt/se'

rule bbtrim:
    input:
        get_fastq
    output:
        fastq = 'bb_trimmed/{sample}-{rep}.fastq.gz'
    params:
        ref = '{}'.format(','.join([os.path.join(snake_dir, i) for i in ['resources/polyA.fa.gz', 'resources/truseq_rna.fa.gz']])),
        options = lambda wildcards: get_program_params(wildcards, program = 'bbtrim', key = 'options')
    log:
        'logs/bbtrim/{sample}-{rep}.log'
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/run_bbtrim.py'
