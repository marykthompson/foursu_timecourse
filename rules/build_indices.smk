rule build_star_index:
    input:
        fasta = 'indices/combo_files/{}_star.fa'.format(config['index_name']),
        gtf = 'indices/combo_files/{}.gtf'.format(config['index_name'])
    output:
        outdir = directory('indices/star_index_{index_name}'.format(index_name = config['index_name']))
    params:
        outdir = 'indices/star_index_{index_name}'.format(index_name = config['index_name']),
        extra='--genomeSAsparseD {} --genomeSAindexNbases {} --limitGenomeGenerateRAM {}'
        .format(config['build_star_params']['genomeSAsparseD'],
        config['build_star_params']['genomeSAindexNbases'],
        config['build_star_params']['limitGenomeGenerateRAM'])
    log:
        'indices/log_files/star_build.log'
    threads: 24
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/build_star_index.py'

rule build_kallisto_index:
    input:
        'indices/combo_files/{}_txts.fa'.format(config['index_name'])
    output:
        'indices/kallisto_index/{}.idx'.format(config['index_name'])
    conda:
        '../envs/main.yaml'
    log:
        'indices/log_files/kallisto_build.log'
    shell:
        'kallisto index -i {output} {input} >& {log}'
