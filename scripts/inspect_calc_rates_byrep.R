#Run INSPEcT on Kallisto transcript quantitations
#Calculate first guess rates
#for this version calculate the rates from single replicates,
#so that the replicate values can be used for the AUC calculations
library('INSPEcT')

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

#input files
nas_exon_file <- snakemake@input[['nas_exon_file']]
nas_intron_file <- snakemake@input[['nas_intron_file']]
tot_exon_file <- snakemake@input[['tot_exon_file']]
tot_intron_file <- snakemake@input[['tot_intron_file']]
exp_des_file <- snakemake@input[['exp_des_file']]

#params
labeling_time <- snakemake@params[['labeling_time']]

#Read in nascent and total exonic and intronic counts and convert to matrix
nasexon_ma <- as.matrix(read.csv(nas_exon_file, row.names = 'gene'))
nasintron_ma <- as.matrix(read.csv(nas_intron_file, row.names = 'gene'))
totexon_ma <- as.matrix(read.csv(tot_exon_file, row.names = 'gene'))
totintron_ma <- as.matrix(read.csv(tot_intron_file, row.names = 'gene'))

#Read in experimental design (tpts, then reps): e.g. t1_1 t2_1, t3_1, t1_2, t2_2, t3_2
#or e.g.: wt_1, syp_1, wt_2, syp_2,...
#I have changed this from timepoints to conditions in the current version
exp_des <- read.csv(exp_des_file)$conditions
num_tpts <- length(unique(exp_des))
tpts <- exp_des[1: num_tpts]
print('exp_des', exp_des)
print('num_tpts', num_tpts)
print('tpts', tpts)

i = 1
r = c(1:num_tpts)
while (i <= length(exp_des)/num_tpts) {
  nasexon_r <- nasexon_ma[, r]
  nasexon_ma_r <- nasexon_ma[,r]
  nasintron_ma_r <- nasintron_ma[,r]
  totexon_ma_r <- totexon_ma[,r]
  totintron_ma_r <- totintron_ma[,r]
  exp_des_r <- exp_des[r]

  nas_L_r <- list('exonsAbundances' = nasexon_ma_r, 'intronsAbundances' = nasintron_ma_r)
  tot_L_r <- list('exonsAbundances' = totexon_ma_r, 'intronsAbundances' = totintron_ma_r)

  nasExp_plgem_r <-quantifyExpressionsFromTrAbundance(trAbundaces = nas_L_r,
                                                      experimentalDesign = exp_des_r)

  totExp_plgem_r <-quantifyExpressionsFromTrAbundance(trAbundaces = tot_L_r,
                                                      experimentalDesign = exp_des_r)

  nascentInspObj_r <-newINSPEcT(tpts = tpts
                                ,labeling_time = labeling_time
                                ,nascentExpressions = nasExp_plgem_r
                                ,matureExpressions = totExp_plgem_r)

  write.csv(ratesFirstGuess(nascentInspObj_r, 'synthesis'), file = paste('inspect/synthesis_', i, '.csv', sep=''))
  write.csv(ratesFirstGuess(nascentInspObj_r, 'degradation'), file = paste('inspect/degradation_', i, '.csv', sep=''))
  write.csv(ratesFirstGuess(nascentInspObj_r, 'processing'), file = paste('inspect/processing_', i, '.csv', sep=''))
  write.csv(ratesFirstGuess(nascentInspObj_r, 'total'), file = paste('inspect/total_', i, '.csv', sep=''))
  write.csv(ratesFirstGuess(nascentInspObj_r, 'preMRNA'), file = paste('inspect/premrna_', i, '.csv', sep=''))
  saveRDS(nascentInspObj_r, file = paste('inspect/inspect_data1_', i, '.rds', sep=''))
  i = i + 1
  r = r + num_tpts
}
