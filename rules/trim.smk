def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ['fq1', 'fq2']].dropna()

rule cutadapt_pe:
    input:
        get_fastq
    output:
        fastq1 = 'cutadapt_trimmed/{sample}-{unit}.1.fastq.gz',
        fastq2 = 'cutadapt_trimmed/{sample}-{unit}.2.fastq.gz',
        qc = 'cutadapt_trimmed/{sample}-{unit}.qc.txt'
    params:
        '-a {} {}'.format(config['trimming']['adapter'], config['params']['cutadapt-pe'])
    log:
        'logs/cutadapt/{sample}-{unit}.log'
    wrapper:
        '0.49.0/bio/cutadapt/pe'

rule cutadapt:
    input:
        get_fastq
    output:
        fastq = 'cutadapt_trimmed/{sample}-{unit}.fastq.gz',
        qc = 'cutadapt_trimmed/{sample}-{unit}.qc.txt'
    params:
        '-a {} {}'.format(config['trimming']['adapter'], config['params']['cutadapt-se'])
    log:
        'logs/cutadapt/{sample}-{unit}.log'
    wrapper:
        '0.49.0/bio/cutadapt/se'

rule bbtrim:
    input:
        get_fastq
    output:
        fastq = 'bb_trimmed/{sample}-{unit}.fastq.gz'
    params:
        extra = 'ref={} {}'.format(','.join([os.path.join(snake_dir, i) for i in config['trimming']['contaminant_files']]), config['params']['bbtrim'])
    log:
        'logs/bbtrim/{sample}-{unit}.log'
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/run_bbtrim.py'
