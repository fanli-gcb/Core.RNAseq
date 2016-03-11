#!/usr/bin/Rscript

# load required packages
library(ggplot2)
library(reshape2)
library(plyr)

args <- commandArgs(T)
if (length(args) < 3) {
	cat("USAGE: ./summ_duplication_stats.R dir out_txt_fn out_pdf_fn \n")
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
	fn <- sprintf("%s/%s.EstimatedLibraryComplexity.txt", dir, sample_id)
	if (!file.exists(fn)) {
		samples <- setdiff(samples, sample_id)
		next
	}
	tmp <- read.table(fn, nrows=1, header=T, quote="", fill=T, blank.lines.skip=T)
	data <- rbind(data, tmp)
}
rownames(data) <- samples
write.table(data, file=out_txt_fn, quote=F, sep="\t", row.names=T, col.names=T)
data$LIBRARY <- rownames(data)

pdf(out_pdf_fn, width=20)
ggplot(data, aes(x=LIBRARY, y=ESTIMATED_LIBRARY_SIZE)) + geom_bar(stat="identity") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
ggplot(data, aes(x=LIBRARY, y=PERCENT_DUPLICATION)) + geom_bar(stat="identity") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
ggplot(data, aes(x=READ_PAIRS_EXAMINED, y=ESTIMATED_LIBRARY_SIZE)) + geom_point() + theme_classic()


dev.off()





