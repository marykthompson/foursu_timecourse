'''
extract_spike_transcripts.py
Extract the transcript seqs for each SIRV or ERCC to use for building a transcript index
'''
import gffutils
from pyfaidx import Fasta

spike_fasta = snakemake.input['spike_fasta']
spike_gtf = snakemake.input['spike_gtf']
spliced_spikes = snakemake.output['spliced_spikes']
db_out = '%s.db' % spike_fasta.split('.')[0]

#Load genes from the fasta file
genome = Fasta(spike_fasta)
db = gffutils.create_db(spike_gtf, db_out, disable_infer_genes = True, disable_infer_transcripts = True, force = True)

txts = db.features_of_type('transcript')

with open(spliced_spikes, 'w') as g:
    for txt in txts:
        txt_id = txt.id
        chrom = txt.chrom
        strand = txt.strand
        if strand == '-':
            exons = db.children(txt, featuretype = 'exon', order_by = 'start', reverse = True)
        else:
            exons = db.children(txt, featuretype = 'exon', order_by = 'start')
        seq = ''
        for exon in exons:
            exon_seq = exon.sequence(genome, use_strand = True)
            seq += exon.sequence(genome, use_strand = True)
        g.write('>%s\n%s\n' % (txt_id, seq))
