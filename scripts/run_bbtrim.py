'''
run_bbtrim.py
Use bbduk to trim adapters and poly(A) from Quant-Seq reads as recommended by Lexogen.
'''

import os
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

sample = [snakemake.input] if isinstance(snakemake.input, str) else snakemake.input
n = len(sample)
assert n == 1 or n == 2, 'input->sample must have 1 (single-end) or 2 (paired-end) elements.'

outprefix = os.path.dirname(snakemake.output[0]) + '/'

shell(
    'bbduk.sh '
    'in={snakemake.input} '
    'out={snakemake.output.fastq} '
    'ref={snakemake.params.ref} '
    '{snakemake.params.options} '
    '{log}')
