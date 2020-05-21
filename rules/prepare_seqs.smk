'''
Download and unzip the files.
Create combined gtf and fasta files from the genome and spike-ins.
'''

rule download_spikes:
    input:
    output:
        'indices/spike_files/{}.fa'.format(config['spike_name']),
        'indices/spike_files/{}.gtf'.format(config['spike_name'])
    params:
        sirv_zip_address = config['sirv_zip_address'],
        ercc_fasta_address = config['ercc_fasta_address'],
        spike_name = config['spike_name'],
        outdir = 'indices/spike_files',
        include_sirvs = config['include_sirvs'],
        include_erccs = config['include_erccs'],
        ercc_subset = config['ercc_subset']
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/prepare_spike_seqs.py'

rule download_genomic_seqs:
    input:
    output:
        fasta = 'indices/genome_files/{}_genome.fa'.format(config['genome_name']),
        gtf = 'indices/genome_files/{}.gtf'.format(config['genome_name']),
        cdna = 'indices/genome_files/{}_cdna.fa'.format(config['genome_name']),
        ncrna = 'indices/genome_files/{}_ncrna.fa'.format(config['genome_name'])
    params:
        fasta_address = config['genome_fasta_address'],
        gtf_address = config['genome_gtf_address'],
        cdna_address = config['cdna_fasta_address'],
        ncrna_address = config['ncrna_fasta_address'],
        outdir = 'indices/genome_files'
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/prepare_ensembl_seqs.py'

rule combine_genome_and_spike:
    input:
        fastas = ['indices/genome_files/{}_genome.fa'.format(config['genome_name']),
        'indices/spike_files/{}.fa'.format(config['spike_name'])],
        gtfs = ['indices/genome_files/{}.gtf'.format(config['genome_name']),
        'indices/spike_files/{}.gtf'.format(config['spike_name'])]
    output:
        combo_fasta = 'indices/combo_files/{}_star.fa'.format(config['index_name']),
        combo_gtf = 'indices/combo_files/{}.gtf'.format(config['index_name'])
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/combine_genomes.py'

rule extract_intron_seqs:
    input:
        combo_fasta = 'indices/combo_files/{}_star.fa'.format(config['index_name']),
        combo_gtf = 'indices/combo_files/{}.gtf'.format(config['index_name'])
    params:
        flanking_nt = 30
    output:
        intron_fasta = 'indices/combo_files/{}_introns.fa'.format(config['index_name']),
        gffutils_db = 'indices/combo_files/{}.db'.format(config['index_name']),
        txt_2_gene_file = 'indices/combo_files/{}_txt2gene.txt'.format(config['index_name'])
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/extract_introns.py'

#Get the spliced versions of the SIRVs to make the kallisto index
rule extract_spike_transcripts:
    input:
        spike_fasta = 'indices/spike_files/{}.fa'.format(config['spike_name']),
        spike_gtf = 'indices/spike_files/{}.gtf'.format(config['spike_name']),
    output:
        spliced_spikes = 'indices/spike_files/{}_spliced.fa'.format(config['spike_name'])
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/extract_spike_transcripts.py'


rule combine_transcripts_and_introns:
    input:
        cdna = 'indices/genome_files/{}_cdna.fa'.format(config['genome_name']),
        ncrna = 'indices/genome_files/{}_ncrna.fa'.format(config['genome_name']),
        spike = 'indices/spike_files/{}_spliced.fa'.format(config['spike_name']),
        intron = 'indices/combo_files/{}_introns.fa'.format(config['index_name'])
    output:
        final_fasta = 'indices/combo_files/{}_txts.fa'.format(config['index_name'])
    shell:
        'cat {input.cdna} {input.ncrna} {input.intron} {input.spike} > {output.final_fasta}'
