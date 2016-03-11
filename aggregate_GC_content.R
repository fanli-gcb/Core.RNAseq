#!/usr/bin/Rscript

# load required packages
library(reshape2)
library(plyr)

args <- commandArgs(T)
if (length(args) < 2) {
	cat("USAGE: ./aggregate_GC_content.R in_file out_file\n")
	q()
}

in_fn <- args[1]
out_fn <- args[2]

data <- read.table(in_fn, header=T, as.is=T, sep="\t", comment.char="")
data$gene_id <- unlist(lapply(data$X9_usercol, function(x) gsub("gene_id ", "", unlist(strsplit(x, ";"))[1])))
data$transcript_id <- unlist(lapply(data$X9_usercol, function(x) gsub(" transcript_id ", "", unlist(strsplit(x, ";"))[2])))

# find longest transcript for each gene
agg <- aggregate(X18_seq_len~transcript_id+gene_id, data, sum)
longest_transcript <- ddply(agg, .(gene_id), function(x) {
	x[which.max(x$X18_seq_len), ]
	})

# compute GC content for each transcript
gc <- ddply(data, .(transcript_id), function(x) {
	sum(x$X13_num_C, x$X14_num_G)/sum(x$X18_seq_len) })
rownames(gc) <- gc$transcript_id

longest_transcript$pct_gc <- gc[longest_transcript$transcript_id, "V1"]
colnames(longest_transcript) <- c("transcript_id", "gene_id", "length", "pct_gc")

write.table(longest_transcript, file=out_fn, quote=F, sep="\t", row.names=F, col.names=T)
