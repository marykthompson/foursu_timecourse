#Rate modeling:
#Feed the Kallisto transcript quant into INSPEcT

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
        excluded_exps = config['exps_excluded_from_inspect']
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/prepare_inspect_input.py'

#run INSPEcT
#for troubleshooting, break into 1 and 2
#1 finds the first guess rates
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
