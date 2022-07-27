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
nas_exp_des_file <- snakemake@input[['nas_exp_des_file']]
tot_exp_des_file <- snakemake@input[['tot_exp_des_file']]

#params
labeling_time <- snakemake@params[['labeling_time']]

#Read in nascent and total exonic and intronic counts and convert to matrix
nasexon_ma <- as.matrix(read.csv(nas_exon_file, row.names = 'gene'))
nasintron_ma <- as.matrix(read.csv(nas_intron_file, row.names = 'gene'))
totexon_ma <- as.matrix(read.csv(tot_exon_file, row.names = 'gene'))
totintron_ma <- as.matrix(read.csv(tot_intron_file, row.names = 'gene'))

#Read in experimental design (tpts, then reps): e.g. t1_1 t2_1, t3_1, t1_2, t2_2, t3_2
#or e.g.: wt_1, syp_1, wt_2, syp_2,...
nas_exp_des <- read.csv(nas_exp_des_file)$conditions
tot_exp_des <- read.csv(tot_exp_des_file)$conditions
num_tpts <- length(unique(nas_exp_des))
tpts <- nas_exp_des[1: num_tpts]

i = 1
r = c(1:num_tpts)
while (i <= length(nas_exp_des)/num_tpts) {
  nasexon_r <- nasexon_ma[, r, drop=FALSE]
  nasexon_ma_r <- nasexon_ma[,r, drop=FALSE]
  nasintron_ma_r <- nasintron_ma[,r, drop=FALSE]
  totexon_ma_r <- totexon_ma[,r, drop=FALSE]
  totintron_ma_r <- totintron_ma[,r, drop=FALSE]
  nas_exp_des_r <- nas_exp_des[r]
  tot_exp_des_r <- tot_exp_des[r]


  nas_L_r <- list('exonsAbundances' = nasexon_ma_r, 'intronsAbundances' = nasintron_ma_r)
  tot_L_r <- list('exonsAbundances' = totexon_ma_r, 'intronsAbundances' = totintron_ma_r)

  nasExp_plgem_r <-quantifyExpressionsFromTrAbundance(trAbundaces = nas_L_r,
                                                      experimentalDesign = nas_exp_des_r)

  totExp_plgem_r <-quantifyExpressionsFromTrAbundance(trAbundaces = tot_L_r,
                                                      experimentalDesign = tot_exp_des_r)

  if (ncol(nasExp_plgem_r$exonsExpressions) != length(unique(nas_exp_des_r))) {
    nas_exonsExpressions2 = t(nasExp_plgem_r$exonsExpressions)
    nas_exonsVariance2 = t(nasExp_plgem_r$exonsVariance)
    nas_intronsExpressions2 = t(nasExp_plgem_r$intronsExpressions)
    nas_intronsVariance2 = t(nasExp_plgem_r$intronsVariance)

    tot_exonsExpressions2 = t(totExp_plgem_r$exonsExpressions)
    tot_exonsVariance2 = t(totExp_plgem_r$exonsVariance)
    tot_intronsExpressions2 = t(totExp_plgem_r$intronsExpressions)
    tot_intronsVariance2 = t(totExp_plgem_r$intronsVariance)

    #add the row names back to the variance matrices
    rownames(nas_exonsVariance2) <- rownames(nas_exonsExpressions2)
    rownames(nas_intronsVariance2) <- rownames(nas_intronsExpressions2)

    rownames(tot_exonsVariance2) <- rownames(tot_exonsExpressions2)
    rownames(tot_intronsVariance2) <- rownames(tot_intronsExpressions2)

    nasExp_plgem_r <- list(exonsExpressions = nas_exonsExpressions2, exonsVariance = nas_exonsVariance2,
                         intronsExpressions = nas_intronsExpressions2, intronsVariance = nas_intronsVariance2)

    totExp_plgem_r <- list(exonsExpressions = tot_exonsExpressions2, exonsVariance = tot_exonsVariance2,
                         intronsExpressions = tot_intronsExpressions2, intronsVariance = tot_intronsVariance2)
  }

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
