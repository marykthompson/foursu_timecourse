def is_single_end(sample, rep):
    return pd.isnull(units.loc[(sample, rep), 'fq2'])

def get_program_params(wildcards, program = '', key = ''):
    '''
    Get the params by libtype. Different libtypes will have different
    mapping strands, etc.
    Program is the software, e.g. kallisto, star
    key is the type of info, e.g. adapter1, se/pe
    '''
    libtype = units.loc[(wildcards.sample, wildcards.rep), 'libtype']
    #pass the key options for the general options specified under se/pe
    if key == 'options':
        if is_single_end(**wildcards):
            key = 'se'
        else:
            key = 'pe'

    extra = config['params'][program][libtype][key]
    return extra
