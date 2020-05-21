'''
run_kallisto.py
Run kallisto quant.

From the snakemake wrappers repo. Author: Joël Simoneau
https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/kallisto/quant.html
'''

from snakemake.shell import shell

# Creating log
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Placeholder for optional parameters
extra = snakemake.params.get('extra', '')

# Allowing for multiple FASTQ files
fastq = snakemake.input.get('fastq')
assert fastq is not None, 'input-> a FASTQ-file is required'
fastq = ' '.join(fastq) if isinstance(fastq, list) else fastq

shell(
    'kallisto quant '
    '{extra} '  # Optional parameters
    '--threads={snakemake.threads} '  # Number of threads
    '--index={snakemake.input.kallisto_index} '  # Input file
    '--output-dir={snakemake.params.outdir} '  # Output directory
    '{fastq} '  # Input FASTQ files
    '{log}'  # Logging
)
