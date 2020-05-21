'''
extract_introns.py
- Extract intronic sequence +30 nt flanking for each intron and create the
gffutils database which will be used for subsequent analyses.
- Also make the transcript to gene file which is needed from some subsequent
programs like DESeq2.
'''

import gffutils
from pyfaidx import Fasta

gtf_file = snakemake.input['combo_gtf']
fasta_file = snakemake.input['combo_fasta']
db_out = snakemake.output['gffutils_db']
intron_fasta = snakemake.output['intron_fasta']
txt_2_gene_file = snakemake.output['txt_2_gene_file']
flanking_nt = snakemake.params['flanking_nt']

def get_txt_info(feat, db):
    gene = next(db.parents(feat, featuretype = 'gene'))
    try:
        symbol = db[gene].attributes['gene_name'][0]
    except KeyError:
        symbol = gene.id
    return gene.id, symbol

db = gffutils.create_db(gtf_file, db_out, disable_infer_genes = True, disable_infer_transcripts = True,
force = True, merge_strategy = 'create_unique')
introns = list(db.create_introns())
#Add introns back to the database
db.update(introns, disable_infer_genes = True, disable_infer_transcripts = True, merge_strategy = 'create_unique')

#Load genes from the fasta file
genome = Fasta(fasta_file)

all_txts = db.features_of_type('transcript')

with open(txt_2_gene_file, 'w') as k:
    with open(intron_fasta, 'w') as g:
        for txt in all_txts:
            txt_id = txt.id
            chrom = txt.chrom
            gene, symbol = get_txt_info(txt, db)
            k.write('%s\t%s\t%s\n' % (txt_id, gene, symbol))

            strand = txt.strand
            if strand == '-':
                introns = [i for i in db.children(txt, featuretype = 'intron', order_by = 'start', reverse = True)]
            else:
                introns = [i for i in db.children(txt, featuretype = 'intron', order_by = 'start')]

            j = 1
            for intron in introns:
                name = '%s_I.%s' % (txt_id, j)
                start = intron.start - flanking_nt - 1
                end = intron.end + flanking_nt
                if strand == '+':
                    seq = genome[chrom][start: end].seq
                else:
                    seq = genome[chrom][start: end].complement.seq[::-1]
                g.write('>%s\n%s\n' % (name, seq))
                k.write('%s\t%s\t%s\n' % (name, gene, symbol))
                j += 1
