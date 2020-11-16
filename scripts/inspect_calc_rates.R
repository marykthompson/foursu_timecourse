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
exp_des_file <- snakemake@input[['exp_des_file']]

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

#Read in nascent and mature exonic and intronic counts and convert to matrix
nasexon_ma <- as.matrix(read.csv(nas_exon_file, row.names = 'gene'))
nasintron_ma <- as.matrix(read.csv(nas_intron_file, row.names = 'gene'))
totexon_ma <- as.matrix(read.csv(tot_exon_file, row.names = 'gene'))
totintron_ma <- as.matrix(read.csv(tot_intron_file, row.names = 'gene'))

#Read in experimental design (tpts, then reps): e.g. t1_1 t2_1, t3_1, t1_2, t2_2, t3_2
#or e.g.: wt_1, syp_1, wt_2, syp_2,...
exp_des <- read.csv(exp_des_file)$conditions
nasL <- list('exonsAbundances' = nasexon_ma, 'intronsAbundances' = nasintron_ma)
totL <- list('exonsAbundances' = totexon_ma, 'intronsAbundances' = totintron_ma)

#Get variances from the replicates
nasExp_plgem<-quantifyExpressionsFromTrAbundance(trAbundaces = nasL,
                                                 experimentalDesign = exp_des)

totExp_plgem<-quantifyExpressionsFromTrAbundance(trAbundaces = totL,
                                                 experimentalDesign = exp_des)

#We will assume from the previous step that the exp_des is exactly each timepoint * n reps
num_tpts <- length(unique(exp_des))
tpts <- exp_des[1: num_tpts]

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

#add the output of the variances here!
saveRDS(nascentInspObj, file = rdata_file)
