'''
prepare_spike_seqs.py
- Download the SIRVs file and fix some of the formatting issues so that it can
be used to build a STAR index.
'''

import subprocess
from zipfile import ZipFile
import os
import shutil
from urllib.request import urlretrieve
from Bio import SeqIO
import gffutils
from Bio.SeqRecord import SeqRecord

def write_gtf(outname, db):
    '''
    Write gtf outfile
    '''
    out_gtf = '{}.gtf'.format(outname)
    with open(out_gtf, 'w') as g:
        for feature in db.all_features():
            g.write(str(feature) + '\n')
    return out_gtf

def download_ERCC_info(address, outdir, new_name, subset = set(), starting_nt = 0, num_Ns = 1000):
    '''
    If there are no SIRVs, start at nt index 0. Else, start at length of SIRV genome
    '''
    outname = os.path.join(outdir, 'erccs')
    download_address = '{}.fasta'.format(outname)
    local_path, headers = urlretrieve(address, download_address)

    ercc_features = []
    ercc_records = SeqIO.parse(download_address, 'fasta')

    this_start = starting_nt
    ercc_combined_seq = ''
    for record in ercc_records:
        if record.id in subset:
            gtf_start = this_start + 1
            gtf_end = this_start + len(record.seq)
            ercc_combined_seq += record.seq + 'N'*num_Ns
            #add an a to the id for the transcript name. gffutils won't infer gene boundaries if gene and transcript id are the same
            gtf_line = '{new_name}\tERCC\texon\t{start}\t{end}\t0\t+\t.\tgene_id "{ID}"; transcript_id "{ID}a";'.format(
            new_name = new_name, ID = record.id, start = gtf_start, end = gtf_end)
            this_start += len(record.seq) + num_Ns
            ercc_features.append(gffutils.feature.feature_from_line(gtf_line))

    ercc_db_file = '{}_db'.format(outname)
    ercc_db = gffutils.create_db(ercc_features, ercc_db_file)
    ercc_gtf = write_gtf(outname, ercc_db)

    os.remove(download_address)
    os.remove(ercc_db_file)

    return ercc_gtf, ercc_combined_seq

def download_SIRV_info(address, outdir, new_name):
    local_path, headers = urlretrieve(address, 'temp.zip')
    print('local_path', local_path)
    print('Zip file downloaded.')
    #Extract files into new dir, move up and delete the default zip name
    extracted_dir = outdir
    with ZipFile('temp.zip', 'r') as zip_file:
        fasta_file = 'SIRV_Set1_Sequences_170612a (ZIP)/SIRVome_isoforms_170612a.fasta'
        gtf_file = 'SIRV_Set1_Sequences_170612a (ZIP)/SIRVome_isoforms_C_170612a.gtf'

        zip_file.extract(fasta_file, extracted_dir)
        zip_file.extract(gtf_file, extracted_dir)

    old_dir, gtf_base = gtf_file.split('/')
    old_dir, fasta_base = fasta_file.split('/')

    new_gtf_path = os.path.join(extracted_dir, gtf_base)
    new_fasta_path = os.path.join(extracted_dir, fasta_base)
    renamed_gtf = os.path.join(extracted_dir, '{}_renamed.gtf'.format(gtf_base.split('.')[0]))

    shutil.move(os.path.join(extracted_dir, fasta_file), new_fasta_path)
    shutil.move(os.path.join(extracted_dir, gtf_file), new_gtf_path)
    shutil.rmtree(os.path.join(extracted_dir, old_dir))

    #replace the chromosome name with the new name
    sirv_features = []
    with open(new_gtf_path, 'r') as f:
        with open(renamed_gtf, 'w') as g:
            for line in f:
                fields = line.split('\t')
                fields[0] = new_name
                g.write('\t'.join(fields))

    print('Done.')

    seq = next(SeqIO.parse(new_fasta_path, "fasta")).seq

    os.remove('temp.zip')
    os.remove(new_fasta_path)
    os.remove(new_gtf_path)

    return renamed_gtf, seq

def main():
    num_Ns = 1000
    new_seq = ''
    gtfs = []
    if snakemake.params['include_sirvs']:
        #return the gtf and the seq (as a Bio.Seq object)
        sirv_gtf, sirv_seq = download_SIRV_info(snakemake.params["sirv_zip_address"],
        snakemake.params["outdir"], snakemake.params["spike_name"])
        gtfs.append(sirv_gtf)
        new_seq += sirv_seq

    else:
        #If not attaching to the SIRV file, start with 5' ends
        new_seq += 'N'*num_Ns

    if snakemake.params['include_erccs']:
        ercc_subset = set(snakemake.params["ercc_subset"])
        ercc_gtf, ercc_seq = download_ERCC_info(snakemake.params["ercc_fasta_address"],
        snakemake.params["outdir"], snakemake.params["spike_name"], subset = ercc_subset,
        starting_nt = len(new_seq), num_Ns = num_Ns)
        gtfs.append(ercc_gtf)
        new_seq += ercc_seq

    #write combined files
    combo_outname  = os.path.join(snakemake.params["outdir"], snakemake.params["spike_name"])

    #write combined gtf
    with open("{}.gtf".format(combo_outname), 'wb') as g:
        for i in gtfs:
            with open(i, 'rb') as this_gtf:
                shutil.copyfileobj(this_gtf, g)
                os.remove(i)

    #write combined fasta
    combo_record = SeqRecord(new_seq, id = snakemake.params["spike_name"], description = "")
    combo_outname  = os.path.join(snakemake.params["outdir"], snakemake.params["spike_name"])
    with open("{}.fa".format(combo_outname), "w") as g:
        SeqIO.write(combo_record, g, "fasta")

if __name__ == '__main__':
    main()
