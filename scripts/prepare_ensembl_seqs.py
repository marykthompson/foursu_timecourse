'''
prepare_ensembl_seqs.py
Download the Ensembl sequences
'''

import subprocess
from zipfile import ZipFile
import os
import shutil
from urllib.request import urlretrieve
import gzip

def gunzip(infile):
    '''Unzip a file with .gz extension. Will remove extension in outfile'''
    with gzip.open(infile, 'rb') as f_in:
        with open(infile.rstrip('.gz'), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(infile)

outdir = snakemake.params["outdir"]
if os.path.exists(outdir):
    shutil.rmtree(outdir)
os.makedirs(outdir, exist_ok = True)

fasta_address = snakemake.params['fasta_address']
gtf_address = snakemake.params['gtf_address']
cdna_address = snakemake.params['cdna_address']
ncrna_address = snakemake.params['ncrna_address']

fasta = "{}.gz".format(snakemake.output["fasta"])
gtf = "{}.gz".format(snakemake.output["gtf"])
cdna = "{}.gz".format(snakemake.output['cdna'])
ncrna = "{}.gz".format(snakemake.output['ncrna'])
#download genomic fasta
#note that ensembl actually uses the unix checksum, not md5
print('downloading fasta')
local_path, headers = urlretrieve(fasta_address, fasta)

#download gtf
print('downloading gtf')
local_path, headers = urlretrieve(gtf_address, gtf)

#download cdna
print('downloading cdna')
local_path, headers = urlretrieve(cdna_address, cdna)

#download ncrna
print('downloading ncrna')
local_path, headers = urlretrieve(ncrna_address, ncrna)

print('Downloaded genome files.')
gunzip(fasta)
gunzip(gtf)
gunzip(cdna)
gunzip(ncrna)

print('Unzipped.')
print('Done.')
