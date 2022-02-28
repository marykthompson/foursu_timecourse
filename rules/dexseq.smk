dexseq_scripts_dir = config['dexseq_scripts_dir']

#Remove the problematic genes with transcripts on two strands
rule clean_gtf:
    input:
        gtf = config['gtf_file'],
        problematic_ids = config['problematic_ids']
    output:
        clean_gtf = 'dexseq/cleaned.gtf'
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/clean_gtf.py'

rule bin_gtf:
    input:
        gtf = 'dexseq/cleaned.gtf'
    output:
        gff = 'dexseq/binned.gff'
    params:
        scriptname = os.path.join(dexseq_scripts_dir, 'dexseq_prepare_annotation.py')
    conda:
        '../envs/main.yaml'
    shell:
        'python {params.scriptname} {input.gtf} {output.gff}'

rule assign_reads:
    input:
        bam = 'star/{sample}-{rep}/Aligned.out.bam',
        gff_file = 'dexseq/binned.gff'
    output:
        exon_counts = 'dexseq/{sample}-{rep}.txt'
    params:
        scriptname = os.path.join(dexseq_scripts_dir, 'dexseq_count.py'),
        extra = lambda wildcards: get_program_params(wildcards, program = 'dexseq', key = 'options')
    conda:
        '../envs/main.yaml'
    shell:
        'python {params.scriptname} {params.extra} {input.gff_file} {input.bam} {output.exon_counts}'
