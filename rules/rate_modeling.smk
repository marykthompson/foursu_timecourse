#Rate modeling:
#Feed the Kallisto transcript quant into INSPEcT

###Prepare a filtered version of the units.tsv file to include experiments used for inspect###
filt_units = units.copy()
filt_units['exp_type'] = filt_units['sample'].apply(lambda x: x.split('_')[0])
filt_units['timepoint'] = filt_units['sample'].apply(lambda x: x.split('_')[1])
filt_units['expname'] = 't' + filt_units['timepoint'] + '_' + filt_units['unit']

#convert replicate number to integer
filt_units['rep_num'] = filt_units['unit'].apply(lambda x: x.split('rep')[1])

#convert these to numeric so that we can sort them
filt_units['rep_num'] = pd.to_numeric(filt_units['rep_num'], errors = 'coerce')
#drop non-numerical replicates
filt_units.dropna(subset = ['rep_num'], inplace = True)
filt_units['timepoint'] = pd.to_numeric(filt_units['timepoint'])

#drop the experiments that shouldn't be included in the analysis
excluded_exps = config['exps_excluded_from_inspect']
filt_units = filt_units[~filt_units['sample'].isin(excluded_exps)].copy()

def get_bam_files(wildcards):
    '''Get the bam files corresponding to the pd and input files'''
    filt_units_input = filt_units[filt_units['sample'].apply(lambda x: x.startswith('input'))].copy()
    filt_units_pd = filt_units[filt_units['sample'].apply(lambda x: x.startswith('pd'))].copy()

    d = {'nascent_files':expand('star/{unit.sample}-{unit.unit}/Aligned.out.bam', unit=filt_units_pd.itertuples()),
    'total_files': expand('star/{unit.sample}-{unit.unit}/Aligned.out.bam', unit=filt_units_input.itertuples())}

    return d

#Prepare INSPEcT input files from the Kallisto quantification.
rule prepare_inspect_input:
    input:
        gene_quant_file = 'results/gene_quantification/summary_abundance_by_gene.csv'
    output:
        exp_des_file = 'inspect/expDes.csv',
        nas_exon_file = 'inspect/nas_exon_tpm.csv',
        nas_intron_file = 'inspect/nas_intron_tpm.csv',
        tot_exon_file = 'inspect/tot_exon_tpm.csv',
        tot_intron_file = 'inspect/tot_intron_tpm.csv'
    params:
        tpts = config['foursu_timepoints'],
        remove_spike_inspect = config['remove_spike_inspect'],
        excluded_exps = config['exps_excluded_from_inspect'],
        primary_col = 'TPM_intron',
        mature_col = 'TPM_exon'
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/prepare_inspect_input.py'

#run INSPEcT
#for troubleshooting, break into 1 and 2
#1 finds the first guess rates
#write the expDes file to input into R

#Try this as an alternate path where data is processed from bam file.
#Put all in folder 'inspect2/'
rule make_expDes_file:
    params:
        tpt_dict = config['foursu_timepoints']
    output:
        'inspect2/expDes2.csv'
    run:
        #get the expDes array to run INSPEcT
        tpt_array = [params.tpt_dict[i] for i in filt_units['timepoint'].unique() if i in params.tpt_dict]
        reps = sorted(filt_units['rep_num'].dropna().unique())
        num_reps = len(reps)
        a = tpt_array * num_reps
        tdf = pd.DataFrame(a, columns=['timepoints'])
        tdf.to_csv(output[0], index=False)

#Try to get these from the bam files directly here.
rule run_inspect1_frombam:
    input:
        unpack(get_bam_files),
        exp_des_file = 'inspect2/expDes2.csv'
    output:
        synth_file = 'inspect2/synth_rates.csv',
        deg_file = 'inspect2/deg_rates.csv',
        proc_file = 'inspect2/proc_rates.csv',
        tot_file = 'inspect2/tot_levels.csv',
        premrna_file = 'inspect2/premrna_levels.csv',
        rdata_file = 'inspect2/inspect_data1.rds'
    params:
        labeling_time = config['labeling_time'],
        txt_db = config['txt_db_R']
    conda:
        '../envs/inspect.yaml'
    log:
        'logs/inpsect/inspect1_bam.log'
    script:
        '../scripts/inspect_calc_rates_bam.R'
    #script:
        #'../scripts/inspect_calc_rates_bam.R'
#2 models the rates and tests the goodness of fit.
rule run_inspect2_frombam:
    input:
        rdata_file = 'inspect2/inspect_data1.rds'
    output:
        synth_mod_file = 'inspect2/synth_model.csv',
        deg_mod_file = 'inspect2/deg_model.csv',
        proc_mod_file = 'inspect2/proc_model.csv',
        tot_mod_file = 'inspect2/tot_model.csv',
        premrna_mod_file = 'inspect2/premrna_model.csv',
        rate_pval_file = 'inspect2/rate_pvals.csv',
        gene_class_file = 'inspect2/gene_class.csv',
        chisq_file = 'inspect2/chisq_fit.csv',
        rdata_file = 'inspect2/inspect_data2.rds'
    params:
        labeling_time = config['labeling_time'],
        model_subset = config['model_subset']
    conda:
        '../envs/inspect.yaml'
    log:
        'logs/inpsect/inspect2_bam.log'
    script:
        '../scripts/inspect_model_rates.R'

#run inspect with non-functional mode (not assuming impulse or sigmoid functions)
rule run_inspect2_frombam_nf:
    input:
        rdata_file = 'inspect2/inspect_data1.rds'
    output:
        synth_mod_file = 'inspect2_nf/synth_model.csv',
        deg_mod_file = 'inspect2_nf/deg_model.csv',
        proc_mod_file = 'inspect2_nf/proc_model.csv',
        tot_mod_file = 'inspect2_nf/tot_model.csv',
        premrna_mod_file = 'inspect2_nf/premrna_model.csv',
        rate_pval_file = 'inspect2_nf/rate_pvals.csv',
        gene_class_file = 'inspect2_nf/gene_class.csv',
        rdata_file = 'inspect2_nf/inspect_data2.rds'
    params:
        labeling_time = config['labeling_time'],
        model_subset = config['model_subset']
    conda:
        '../envs/inspect.yaml'
    log:
        'logs/inpsect/inspect2_bam_nf.log'
    script:
        '../scripts/inspect_model_rates_NF.R'

rule run_inspect2_nf:
    input:
        rdata_file = 'inspect/inspect_data1.rds'
    output:
        synth_mod_file = 'inspect_nf/synth_model.csv',
        deg_mod_file = 'inspect_nf/deg_model.csv',
        proc_mod_file = 'inspect_nf/proc_model.csv',
        tot_mod_file = 'inspect_nf/tot_model.csv',
        premrna_mod_file = 'inspect_nf/premrna_model.csv',
        rate_pval_file = 'inspect_nf/rate_pvals.csv',
        gene_class_file = 'inspect_nf/gene_class.csv',
        rdata_file = 'inspect_nf/inspect_data2.rds'
    params:
        labeling_time = config['labeling_time'],
        model_subset = config['model_subset']
    conda:
        '../envs/inspect.yaml'
    log:
        'logs/inpsect/inspect2_nf.log'
    script:
        '../scripts/inspect_model_rates_NF.R'

rule run_inspect1:
    input:
        exp_des_file = 'inspect/expDes.csv',
        nas_exon_file = 'inspect/nas_exon_tpm.csv',
        nas_intron_file = 'inspect/nas_intron_tpm.csv',
        tot_exon_file = 'inspect/tot_exon_tpm.csv',
        tot_intron_file = 'inspect/tot_intron_tpm.csv'
    output:
        synth_file = 'inspect/synth_rates.csv',
        deg_file = 'inspect/deg_rates.csv',
        proc_file = 'inspect/proc_rates.csv',
        tot_file = 'inspect/tot_levels.csv',
        premrna_file = 'inspect/premrna_levels.csv',
        rdata_file = 'inspect/inspect_data1.rds'
    params:
        labeling_time = config['labeling_time'],
    conda:
        '../envs/inspect.yaml'
    log:
        'logs/inpsect/inspect1.log'
    script:
        '../scripts/inspect_calc_rates.R'

#2 models the rates and tests the goodness of fit.
rule run_inspect2:
    input:
        rdata_file = 'inspect/inspect_data1.rds'
    output:
        synth_mod_file = 'inspect/synth_model.csv',
        deg_mod_file = 'inspect/deg_model.csv',
        proc_mod_file = 'inspect/proc_model.csv',
        tot_mod_file = 'inspect/tot_model.csv',
        premrna_mod_file = 'inspect/premrna_model.csv',
        rate_pval_file = 'inspect/rate_pvals.csv',
        gene_class_file = 'inspect/gene_class.csv',
        chisq_file = 'inspect/chisq_fit.csv',
        rdata_file = 'inspect/inspect_data2.rds'
    params:
        labeling_time = config['labeling_time'],
        model_subset = config['model_subset']
    conda:
        '../envs/inspect.yaml'
    log:
        'logs/inpsect/inspect2.log'
    script:
        '../scripts/inspect_model_rates.R'
