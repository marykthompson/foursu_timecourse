#Run INSPEcT on Kallisto transcript quantitations
#Model rates and calculate p-values for rate change and chisq fit
library('INSPEcT')

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

rdata_file_in <- snakemake@input[['rdata_file']]

#output files for modeled rates
synth_mod_file <- snakemake@output[['synth_mod_file']]
deg_mod_file <- snakemake@output[['deg_mod_file']]
proc_mod_file <- snakemake@output[['proc_mod_file']]
tot_mod_file <- snakemake@output[['tot_mod_file']]
premrna_mod_file <- snakemake@output[['premrna_mod_file']]

#output file for the rate p-values
rate_pval_file <- snakemake@output[['rate_pval_file']]
gene_class_file <- snakemake@output[['gene_class_file']]
chisq_file <- snakemake@output[['chisq_file']]
rdata_file_out <- snakemake@output[['rdata_file']]

#params
labeling_time <- snakemake@params[['labeling_time']]
model_subset <- snakemake@params[['model_subset']]
model_rates_nonfunctional <- snakemake@params[['model_rates_nonfunctional']]
model_seed <- snakemake@params[['model_seed']]

print(paste0('model_subset', model_subset))
if (model_subset == TRUE) {
  print('model subset true!')
}

nascentInspObj <- readRDS(rdata_file_in)

#model rates for either a subset or the whole dataset
if (model_subset == TRUE) {
  nascentInspObj_sub <- nascentInspObj[1:20]
} else {
  nascentInspObj_sub <- nascentInspObj[1:nGenes(nascentInspObj)]
}

print('modeling rates')
#model rates using either functional or non-functional form
#modelRatesNF doesn't have the seed option
if (model_rates_nonfunctional == TRUE) {
  nascentInspObj_sub <- modelRatesNF(nascentInspObj_sub)
} else {
  nascentInspObj_sub <- modelRates(nascentInspObj_sub, seed=model_seed)
}

print('writing modeled rates')
write.csv(viewModelRates(nascentInspObj_sub, 'synthesis'), file = synth_mod_file)
write.csv(viewModelRates(nascentInspObj_sub, 'processing'), file = proc_mod_file)
write.csv(viewModelRates(nascentInspObj_sub, 'degradation'), file = deg_mod_file)
write.csv(viewModelRates(nascentInspObj_sub, 'total'), file = tot_mod_file)
write.csv(viewModelRates(nascentInspObj_sub, 'preMRNA'), file = premrna_mod_file)

print('writing supplementary files')
write.csv(ratePvals(nascentInspObj_sub), file = rate_pval_file)
write.csv(geneClass(nascentInspObj_sub), file = gene_class_file)
write.csv(chisqmodel(nascentInspObj_sub), file = chisq_file)

saveRDS(nascentInspObj_sub, file = rdata_file_out)
