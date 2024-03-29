#Rate modeling:
#Feed the Kallisto transcript quant into INSPEcT
import sys

#Prepare INSPEcT input files from the Kallisto quantification.
rule prepare_inspect_input:
    input:
        gene_quant_file = 'results/gene_quantification/summary_abundance_by_gene_filtered.csv'
    output:
        nas_exp_des_file = 'inspect/nas_expDes.csv',
        nas_exon_file = 'inspect/nas_exon_tpm.csv',
        nas_intron_file = 'inspect/nas_intron_tpm.csv',
        tot_exp_des_file = 'inspect/tot_expDes.csv',
        tot_exon_file = 'inspect/tot_exon_tpm.csv',
        tot_intron_file = 'inspect/tot_intron_tpm.csv'
    params:
        condition_mapping = config['foursu_condition_mapping'],
        excluded_exps = config['exps_excluded_from_inspect'],
        analysis_type = config['kinetic_analysis_type'],
        intronic_col = 'intronic_tpm_recalc',
        exonic_col = 'exonic_tpm_recalc'
    conda:
        '../envs/main.yaml'
    script:
        '../scripts/prepare_inspect_input.py'

#run INSPEcT
#break into 1 and 2
#1 finds the first guess rates
#& writes the expDes file to input into R
rule run_inspect1:
    input:
        nas_exp_des_file = 'inspect/nas_expDes.csv',
        tot_exp_des_file = 'inspect/tot_expDes.csv',
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
        synth_var_file = 'inspect/synth_var.csv',
        total_var_file = 'inspect/total_var.csv',
        premrna_var_file = 'inspect/premrna_var.csv',
        rdata_file = 'inspect/inspect_data1.rds'
    params:
        labeling_time = config['labeling_time'],
    conda:
        '../envs/inspect.yaml'
    log:
        'logs/inspect/inspect1.log'
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
        model_subset = config['model_subset'],
        model_rates_nonfunctional = config['model_rates_nonfunctional'],
        model_seed = config['model_seed']
    conda:
        '../envs/inspect.yaml'
    log:
        'logs/inspect/inspect2.log'
    script:
        '../scripts/inspect_model_rates.R'

#Run INSPEcT with each rep separately
rule run_inspect_byrep:
    input:
        nas_exp_des_file = 'inspect/nas_expDes.csv',
        tot_exp_des_file = 'inspect/tot_expDes.csv',
        nas_exon_file = 'inspect/nas_exon_tpm.csv',
        nas_intron_file = 'inspect/nas_intron_tpm.csv',
        tot_exon_file = 'inspect/tot_exon_tpm.csv',
        tot_intron_file = 'inspect/tot_intron_tpm.csv'
    output:
        list(set(expand('inspect/synthesis_{unit.replicate}.csv', unit = units.itertuples()))),
        list(set(expand('inspect/degradation_{unit.replicate}.csv', unit = units.itertuples()))),
        list(set(expand('inspect/processing_{unit.replicate}.csv', unit = units.itertuples()))),
        list(set(expand('inspect/total_{unit.replicate}.csv', unit = units.itertuples()))),
        list(set(expand('inspect/premrna_{unit.replicate}.csv', unit = units.itertuples()))),
        list(set(expand('inspect/inspect_data1_{unit.replicate}.rds', unit = units.itertuples())))
    params:
        labeling_time = config['labeling_time'],
    conda:
        '../envs/inspect.yaml'
    log:
        'logs/inspect/inspect_byrep.log'
    script:
        '../scripts/inspect_calc_rates_byrep.R'
