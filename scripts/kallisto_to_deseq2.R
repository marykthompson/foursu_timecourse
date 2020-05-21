#Use the tximport library and DESeq2 to do gene-level diffexp analysis
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("tximport")
library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

units_file <- snakemake@params[["units_file"]]
samples_file <- snakemake@params[["samples_file"]]

units <- read.table(units_file, sep = '\t', header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
samples <- read.table(samples_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

#add the condition data to the units df
units <- merge(units, samples, by = "sample")
#add the name of the run, this will match the name in the Kallisto dir
units$run <- paste(units$sample, units$unit, sep = "-")
#add the replicate number
units$rep <- sapply(units$unit, function(x) unlist(strsplit(x, "rep"))[2])

this_contrast = snakemake@params[["contrast"]]

txt_2_gene_file = snakemake@input[["txt_2_gene_file"]]

t2g <- read.table(txt_2_gene_file, sep = "\t", header = FALSE,
col.names =  c("target_id", "gene", "symbol"), stringsAsFactors = FALSE)
t2g$symbol <- NULL

files <- snakemake@input[["infiles"]]
names(files) <- sapply(files, function(x) unlist(strsplit(x, "/"))[2])

#Subset the data by experment:
exp_df <- units[units$condition %in% this_contrast,]
#only run files that are present
exp_df <- exp_df[exp_df$run %in% names(files)]

condition_df <- exp_df[, c("run", "condition", "rep")]
txi <- tximport(files, type="kallisto", tx2gene=t2g)
ddsTxi <- DESeqDataSetFromTximport(txi, colData = condition_df, design = ~ condition)
#https://rdrr.io/bioc/DESeq2/man/DESeq.html
dds <- DESeq(ddsTxi, sfType = snakemake@params[["sf_method"]])
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file = snakemake@output[["table"]])

# store results
svg(snakemake@output[["ma_plot"]])
plotMA(res, ylim=c(-2,2))
dev.off()
