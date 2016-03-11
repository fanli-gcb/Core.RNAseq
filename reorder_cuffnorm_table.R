#!/usr/bin/Rscript

# load required packages

args <- commandArgs(T)
if (length(args) < 3) {
	cat("USAGE: ./reorder_cuffnorm_table.R cuffnorm_table_fn sample_table_fn out_fn\n")
	q()
}

cuffnorm_table_fn <- args[1]
sample_table_fn <- args[2]
out_fn <- args[3]

samples <- read.table(sample_table_fn, header=T, sep="\t", as.is=T, row.names=1)
data <- read.table(cuffnorm_table_fn, header=T, sep="\t")
samples$correct_id <- basename(gsub("/abundances.cxb", "", samples$file))
new <- samples[colnames(data)[-1], "correct_id"]
colnames(data) <- c("tracking_id", new)

write.table(data, file=out_fn, sep="\t", row.names=F, col.names=T, quote=F)
