'''
build_star_index.py
Wrapper to build star index
'''

import os
from snakemake.shell import shell
import shutil

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
outdir = snakemake.output["outdir"] + '/'

if os.path.exists(outdir):
    shutil.rmtree(outdir)
os.mkdir(outdir)

shell(
    "STAR "
    "{snakemake.params.extra} "
    "--runThreadN {snakemake.threads} "
    "--runMode genomeGenerate "
    "--genomeDir {outdir} "
    "--genomeFastaFiles {snakemake.input.fasta} "
    "--sjdbGTFfile {snakemake.input.gtf} "
    "--outStd Log "
    "{log}")

#Move the log file produced by the star build to the star build directory
shutil.move('Log.out', os.path.join(outdir, 'Log.out'))
