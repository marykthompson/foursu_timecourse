#Run INSPEcT on Kallisto transcript quantitations
#Calculate first guess rates
library('INSPEcT')

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

#input files
nas_exon_file <- snakemake@input[['nas_exon_file']]
nas_intron_file <- snakemake@input[['nas_intron_file']]
tot_exon_file <- snakemake@input[['tot_exon_file']]
tot_intron_file <- snakemake@input[['tot_intron_file']]
nas_exp_des_file <- snakemake@input[['nas_exp_des_file']]
tot_exp_des_file <- snakemake@input[['tot_exp_des_file']]

#params
labeling_time <- snakemake@params[['labeling_time']]

#output files for first guess rates
synth_file <- snakemake@output[['synth_file']]
deg_file <- snakemake@output[['deg_file']]
proc_file <- snakemake@output[['proc_file']]
tot_file <- snakemake@output[['tot_file']]
premrna_file <- snakemake@output[['premrna_file']]
synth_var_file <- snakemake@output[['synth_var_file']]
total_var_file <- snakemake@output[['total_var_file']]
premrna_var_file <- snakemake@output[['premrna_var_file']]

rdata_file <- snakemake@output[['rdata_file']]

#Read in nascent and total exonic and intronic counts and convert to matrix
nasexon_ma <- as.matrix(read.csv(nas_exon_file, row.names = 'gene'))
nasintron_ma <- as.matrix(read.csv(nas_intron_file, row.names = 'gene'))
totexon_ma <- as.matrix(read.csv(tot_exon_file, row.names = 'gene'))
totintron_ma <- as.matrix(read.csv(tot_intron_file, row.names = 'gene'))

#Read in experimental design (tpts, then reps): e.g. t1_1 t2_1, t3_1, t1_2, t2_2, t3_2
#or e.g.: wt_1, syp_1, wt_2, syp_2,...
nas_exp_des <- read.csv(nas_exp_des_file)$conditions
tot_exp_des <- read.csv(tot_exp_des_file)$conditions

nasL <- list('exonsAbundances' = nasexon_ma, 'intronsAbundances' = nasintron_ma)
totL <- list('exonsAbundances' = totexon_ma, 'intronsAbundances' = totintron_ma)

#Get variances from the replicates
nasExp_plgem<-quantifyExpressionsFromTrAbundance(trAbundaces = nasL,
                                                 experimentalDesign = nas_exp_des)

totExp_plgem<-quantifyExpressionsFromTrAbundance(trAbundaces = totL,
                                                 experimentalDesign = tot_exp_des)

# There is a bug in INSPEcT in which it reduces one condition matrix to vector.
# Transpose the row/columns in order to fix the matrix subsetting problem which
# becomes a problem with one condition.

if (ncol(nasExp_plgem$exonsExpressions) != length(unique(nas_exp_des))) {
  nas_exonsExpressions2 = t(nasExp_plgem$exonsExpressions)
  nas_exonsVariance2 = t(nasExp_plgem$exonsVariance)
  nas_intronsExpressions2 = t(nasExp_plgem$intronsExpressions)
  nas_intronsVariance2 = t(nasExp_plgem$intronsVariance)

  tot_exonsExpressions2 = t(totExp_plgem$exonsExpressions)
  tot_exonsVariance2 = t(totExp_plgem$exonsVariance)
  tot_intronsExpressions2 = t(totExp_plgem$intronsExpressions)
  tot_intronsVariance2 = t(totExp_plgem$intronsVariance)

  #add the row names back to the variance matrices
  rownames(nas_exonsVariance2) <- rownames(nas_exonsExpressions2)
  rownames(nas_intronsVariance2) <- rownames(nas_intronsExpressions2)

  rownames(tot_exonsVariance2) <- rownames(tot_exonsExpressions2)
  rownames(tot_intronsVariance2) <- rownames(tot_intronsExpressions2)

  nasExp_plgem <- list(exonsExpressions = nas_exonsExpressions2, exonsVariance = nas_exonsVariance2,
                       intronsExpressions = nas_intronsExpressions2, intronsVariance = nas_intronsVariance2)

  totExp_plgem <- list(exonsExpressions = tot_exonsExpressions2, exonsVariance = tot_exonsVariance2,
                       intronsExpressions = tot_intronsExpressions2, intronsVariance = tot_intronsVariance2)
}

#This is only true for timecourses maybe?
#Will this work for multiple steadystate conditions?

num_tpts <- length(unique(nas_exp_des))
tpts <- nas_exp_des[1: num_tpts]
#for steady-state data, timepoints will actually be conditions
nascentInspObj<-newINSPEcT(tpts = tpts
                           ,labeling_time = labeling_time
                           ,nascentExpressions = nasExp_plgem
                           ,matureExpressions = totExp_plgem)

print('writing first guess rates')
#write the rates first guess:
write.csv(ratesFirstGuess(nascentInspObj, 'synthesis'), file = synth_file)
write.csv(ratesFirstGuess(nascentInspObj, 'degradation'), file = deg_file)
write.csv(ratesFirstGuess(nascentInspObj, 'processing'), file = proc_file)
write.csv(ratesFirstGuess(nascentInspObj, 'total'), file = tot_file)
write.csv(ratesFirstGuess(nascentInspObj, 'preMRNA'), file = premrna_file)
write.csv(ratesFirstGuessVar(nascentInspObj, 'synthesis'), file = synth_var_file)
write.csv(ratesFirstGuessVar(nascentInspObj, 'total'), file = total_var_file)
write.csv(ratesFirstGuessVar(nascentInspObj, 'preMRNA'), file = premrna_var_file)

saveRDS(nascentInspObj, file = rdata_file)
