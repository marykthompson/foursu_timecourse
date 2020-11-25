'''
run_star.py
Run STAR alignment.

From the snakemake wrappers repo. Author: Johannes Köster
https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/star/align.html
Modified to used the gunzip -c command to read in gzipped files.
'''

import os
from snakemake.shell import shell

extra = snakemake.params.get('extra', '')
log = snakemake.log_fmt_shell(stdout=True, stderr=True)


sample = [snakemake.input.sample] if isinstance(snakemake.input.sample, str) else snakemake.input.sample
n = len(sample)
assert n == 1 or n == 2, 'input->sample must have 1 (single-end) or 2 (paired-end) elements.'

if sample[0].endswith('.gz'):
    readcmd = '--readFilesCommand gunzip -c'
else:
    readcmd = ''


outprefix = os.path.dirname(snakemake.output[0]) + '/'


shell(
    'STAR '
    '--quantMode GeneCounts --sjdbGTFfile '
    '{snakemake.params.gtf} {snakemake.params.options} '
    '--runThreadN {snakemake.threads} '
    '--genomeDir {snakemake.params.star_index} '
    '--readFilesIn {snakemake.input.sample} '
    '{readcmd} '
    '--outSAMtype BAM Unsorted '
    '--outFileNamePrefix {outprefix} '
    '--outStd Log '
    '{log}')
