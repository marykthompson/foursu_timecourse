#Remove problematic gene IDs from the gtf file so that DEXSeq binning can occur
problem_file = snakemake.input['problematic_ids']
gtf_file = snakemake.input['gtf']
gtf_file_clean = snakemake.output['clean_gtf']

with open(problem_file, 'r') as f:
    genes = [i.strip('\n') for i in f.readlines()]
genes = set(genes)

with open(gtf_file, 'r') as f:
    with open(gtf_file_clean, 'w') as g:
        for line in f.readlines():
            #print(line)
            fields = line.split('\t')
            if len(fields) == 1:
                continue
            gene_id = fields[8].split(';')[0].split(' ')[1].strip('"')
            if gene_id in genes:
                continue
            else:
                g.write(line)
