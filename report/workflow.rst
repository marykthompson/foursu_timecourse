This workflow performs differential expression analysis on single- or paired-end RNA-seq data.
Adapters were removed with `Cutadapt <http://cutadapt.readthedocs.io>`_ or with
`BBDuk <https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/>`_
for Quant-Seq data. Transcript abundance was estimated with
`Kallisto <https://pachterlab.github.io/kallisto/manual>`_. Transcript estimates
were combined to gene-level estimates and RNA synthesis, processing, and decay rates
were estimated with `INSPEcT <https://bioconductor.org/packages/release/bioc/html/INSPEcT.html>`_.
