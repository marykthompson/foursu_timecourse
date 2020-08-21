#Run INSPEcT on Kallisto transcript quantitations
#Calculate first guess rates
library('INSPEcT')
library('GenomicFeatures')

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

#input files
nas_bams <- snakemake@input[['nascent_files']]
total_bams <- snakemake@input[['total_files']]
exp_des_file <- snakemake@input[['exp_des_file']]

#params
labeling_time <- snakemake@params[['labeling_time']]
txt_db_file <- snakemake@params[['txt_db']]

#output files for first guess rates
synth_file <- snakemake@output[['synth_file']]
deg_file <- snakemake@output[['deg_file']]
proc_file <- snakemake@output[['proc_file']]
tot_file <- snakemake@output[['tot_file']]
premrna_file <- snakemake@output[['premrna_file']]

rdata_file <- snakemake@output[['rdata_file']]

#Read in experimental design (tpts, then reps): e.g. t1_1 t2_1, t3_1, t1_2, t2_2, t3_2
exp_des <- read.csv(exp_des_file)$timepoints

#Load the txt_db
txt_db <- loadDb(txt_db_file)

#Quantify expression from the Bam files
nasExp <- quantifyExpressionsFromBAMs(txdb=txt_db
           , BAMfiles=nas_bams
           , by = 'gene'
           , allowMultiOverlap = FALSE
           , strandSpecific = 1
           , experimentalDesign=exp_des)

totExp <- quantifyExpressionsFromBAMs(txdb=txt_db
          , BAMfiles=total_bams
          , by = 'gene'
          , allowMultiOverlap = FALSE
          , strandSpecific = 1
          , experimentalDesign=exp_des)

#We will assume from the previous step that the exp_des is exactly each timepoint * n reps
num_tpts <- length(unique(exp_des))
tpts <- exp_des[1: num_tpts]

#create INSPEcT object
nascentInspObj<-newINSPEcT(tpts = tpts
                           ,labeling_time = labeling_time
                           ,nascentExpressions = nasExp
                           ,matureExpressions = totExp)

print('writing first guess rates')
#write the rates first guess:
write.csv(ratesFirstGuess(nascentInspObj, 'synthesis'), file = synth_file)
write.csv(ratesFirstGuess(nascentInspObj, 'degradation'), file = deg_file)
write.csv(ratesFirstGuess(nascentInspObj, 'processing'), file = proc_file)
write.csv(ratesFirstGuess(nascentInspObj, 'total'), file = tot_file)
write.csv(ratesFirstGuess(nascentInspObj, 'preMRNA'), file = premrna_file)

saveRDS(nascentInspObj, file = rdata_file)
