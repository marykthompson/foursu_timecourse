#Genome plotting:
#Make bigwig files from bam files to use for interactive browsing
#more to follow...

#Sort and index bam files
rule sort_and_index:
    input:
        'star/{sample}-{unit}/Aligned.out.bam'
    output:
        'star/{sample}-{unit}/sorted.bam'
    conda:
        '../envs/main.yaml'
    shell:
        '''
        samtools sort {input} -o {output}
        samtools index {output}
        '''

#Make the bigwig file. Use BPM scaling and make separate files for the plus and minus strands.
#Note that the --filterRNAstrand option assumes the dUTP method,
#so for Ovation libaries, .p.bw files will actually be the minus strand.
#https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
rule make_bigwig:
    input:
        'star/{sample}-{unit}/sorted.bam'
    output:
        plus_strand_bw = 'bigwig/{sample}-{unit}.p.bw',
        minus_strand_bw = 'bigwig/{sample}-{unit}.m.bw'
    conda:
        '../envs/main.yaml'
    shell:
        '''
        bamCoverage -b {input} --filterRNAstrand forward --normalizeUsing BPM -o {output.plus_strand_bw}
        bamCoverage -b {input} --filterRNAstrand reverse --normalizeUsing BPM -o {output.minus_strand_bw}
        '''
