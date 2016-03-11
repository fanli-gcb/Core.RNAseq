#!/usr/bin/Rscript

# load required packages
library(ggplot2)
library(reshape2)
library(plyr)

args <- commandArgs(T)
if (length(args) < 3) {
	cat("USAGE: ./summ_rnaseq_stats.R dir out_txt_fn out_pdf_fn \n")
	q()
}

dir <- args[1]
out_txt_fn <- args[2]
out_pdf_fn <- args[3]

samples <- gsub(".unaligned.bam", "", list.files(sprintf("%s/", dir), pattern=".unaligned.bam"))
trim.whitespace <- function(x) {
	sub("\\s+$", "", sub("^\\s+", "", x))
}


data <- {}
for (sample_id in samples) {
	fn <- sprintf("%s/%s.RnaSeqMetrics.txt", dir, sample_id)
	if (!file.exists(fn)) {
		samples <- setdiff(samples, sample_id)
		next
	}
	tmp <- read.table(fn, nrows=1, header=T, quote="", fill=T, blank.lines.skip=T)
	data <- rbind(data, tmp)
}
rownames(data) <- samples
final <- data[, c("PCT_RIBOSOMAL_BASES", "PCT_CODING_BASES", "PCT_UTR_BASES", "PCT_INTRONIC_BASES", "PCT_INTERGENIC_BASES", "PCT_MRNA_BASES")]
colnames(final) <- c("rRNA", "CDS", "UTR", "intron", "intergenic", "exon")
write.table(final, file=out_txt_fn, quote=F, sep="\t", row.names=T, col.names=T)


final2 <- final[, c(1,4,5,6)]
final2$ID <- rownames(final2)
df <- melt(final2, id.vars="ID")
#df$variable <- factor(df$variable, levels=c("exon", "intron", "rRNA", "intergenic"))
df$variable <- factor(df$variable, levels=c("intergenic", "rRNA", "intron", "exon"))

pdf(out_pdf_fn, width=20)
ggplot(df, aes(x=ID, y=value, fill=variable)) + geom_bar(stat="identity", aes(order=desc(variable))) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_brewer(palette="Set1")
dev.off()


