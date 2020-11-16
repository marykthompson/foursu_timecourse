## RSEQC

rule rseqc_gtf2bed:
    input:
        config['gtf_file']
    output:
        bed = 'qc/rseqc/annotation.bed',
        db = temp('qc/rseqc/annotation.db')
    log:
        'logs/rseqc_gtf2bed.log'
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/gtf2bed.py'


rule rseqc_junction_annotation:
    input:
        bam = 'star/{sample}-{rep}/Aligned.out.bam',
        bed = 'qc/rseqc/annotation.bed'
    output:
        'qc/rseqc/{sample}-{rep}.junctionanno.junction.bed'
    priority: 1
    log:
        'logs/rseqc/rseqc_junction_annotation/{sample}-{rep}.log'
    params:
        extra = r'-q 255',  # STAR uses 255 as a score for unique mappers
        prefix= 'qc/rseqc/{sample}-{rep}.junctionanno'
    conda:
        '../envs/main.yaml'
    shell:
        'junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} '
        '> {log[0]} 2>&1'


rule rseqc_junction_saturation:
    input:
        bam = 'star/{sample}-{rep}/Aligned.out.bam',
        bed = 'qc/rseqc/annotation.bed'
    output:
        'qc/rseqc/{sample}-{rep}.junctionsat.junctionSaturation_plot.pdf'
    priority: 1
    log:
        'logs/rseqc/rseqc_junction_saturation/{sample}-{rep}.log'
    params:
        extra = r'-q 255',
        prefix = 'qc/rseqc/{sample}-{rep}.junctionsat'
    conda:
        '../envs/main.yaml'
    shell:
        'junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} '
        '> {log} 2>&1'


rule rseqc_stat:
    input:
        'star/{sample}-{rep}/Aligned.out.bam',
    output:
        'qc/rseqc/{sample}-{rep}.stats.txt'
    priority: 1
    log:
        'logs/rseqc/rseqc_stat/{sample}-{rep}.log'
    conda:
        '../envs/main.yaml'
    shell:
        'bam_stat.py -i {input} > {output} 2> {log}'


rule rseqc_infer:
    input:
        bam = 'star/{sample}-{rep}/Aligned.out.bam',
        bed = 'qc/rseqc/annotation.bed'
    output:
        'qc/rseqc/{sample}-{rep}.infer_experiment.txt'
    priority: 1
    log:
        'logs/rseqc/rseqc_infer/{sample}-{rep}.log'
    conda:
        '../envs/main.yaml'
    shell:
        'infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}'


rule rseqc_innerdis:
    input:
        bam = 'star/{sample}-{rep}/Aligned.out.bam',
        bed = 'qc/rseqc/annotation.bed'
    output:
        'qc/rseqc/{sample}-{rep}.inner_distance_freq.inner_distance.txt'
    priority: 1
    log:
        'logs/rseqc/rseqc_innerdis/{sample}-{rep}.log'
    params:
        prefix = 'qc/rseqc/{sample}-{rep}.inner_distance_freq'
    conda:
        '../envs/main.yaml'
    shell:
        'inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1'


rule rseqc_readdis:
    input:
        bam = 'star/{sample}-{rep}/Aligned.out.bam',
        bed = 'qc/rseqc/annotation.bed'
    output:
        'qc/rseqc/{sample}-{rep}.readdistribution.txt'
    priority: 1
    log:
        'logs/rseqc/rseqc_readdis/{sample}-{rep}.log'
    conda:
        '../envs/main.yaml'
    shell:
        'read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}'


rule rseqc_readdup:
    input:
        'star/{sample}-{rep}/Aligned.out.bam'
    output:
        'qc/rseqc/{sample}-{rep}.readdup.DupRate_plot.pdf'
    priority: 1
    log:
        'logs/rseqc/rseqc_readdup/{sample}-{rep}.log'
    params:
        prefix = 'qc/rseqc/{sample}-{rep}.readdup'
    conda:
        '../envs/main.yaml'
    shell:
        'read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1'


rule rseqc_readgc:
    input:
        'star/{sample}-{rep}/Aligned.out.bam'
    output:
        'qc/rseqc/{sample}-{rep}.readgc.GC_plot.pdf'
    priority: 1
    log:
        'logs/rseqc/rseqc_readgc/{sample}-{rep}.log'
    params:
        prefix = 'qc/rseqc/{sample}-{rep}.readgc'
    conda:
        '../envs/main.yaml'
    shell:
        'read_GC.py -i {input} -o {params.prefix} > {log} 2>&1'


rule multiqc:
    input:
        expand('star/{unit.sample}-{unit.replicate}/Aligned.out.bam', unit=units.itertuples()),
        expand('qc/rseqc/{unit.sample}-{unit.replicate}.junctionanno.junction.bed', unit=units.itertuples()),
        expand('qc/rseqc/{unit.sample}-{unit.replicate}.junctionsat.junctionSaturation_plot.pdf', unit=units.itertuples()),
        expand('qc/rseqc/{unit.sample}-{unit.replicate}.infer_experiment.txt', unit=units.itertuples()),
        expand('qc/rseqc/{unit.sample}-{unit.replicate}.stats.txt', unit=units.itertuples()),
        expand('qc/rseqc/{unit.sample}-{unit.replicate}.inner_distance_freq.inner_distance.txt', unit=units.itertuples()),
        expand('qc/rseqc/{unit.sample}-{unit.replicate}.readdistribution.txt', unit=units.itertuples()),
        expand('qc/rseqc/{unit.sample}-{unit.replicate}.readdup.DupRate_plot.pdf', unit=units.itertuples()),
        expand('qc/rseqc/{unit.sample}-{unit.replicate}.readgc.GC_plot.pdf', unit=units.itertuples()),
        expand('logs/rseqc/rseqc_junction_annotation/{unit.sample}-{unit.replicate}.log', unit=units.itertuples())
    output:
        'qc/multiqc_report.html'
    log:
        'logs/multiqc.log'
    wrapper:
        '0.49.0/bio/multiqc'
