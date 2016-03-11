#!/usr/bin/Rscript

# load required packages
library(RColorBrewer)
library(ggplot2)
library(reshape2)

# prefix and full_file can be comma-separated lists
args <- commandArgs(T)
if (length(args) < 3) {
	cat("USAGE: ./saturation_analysis.R prefix full_file out_pdf\n")
	q()
}

prefixes <- unlist(strsplit(args[1], ","))
full_files <- unlist(strsplit(args[2], ","))
out_pdf <- args[3]

# files are: /Lab_Share/fanli/data/20140401/20140401-1A1.sub_$pct.cufflinks_output/genes.fpkm_tracking

if (length(prefixes) != length(full_files)) {
	cat("Number of prefixes and full_files doesn't match!\n")
	q()
}
pdf(out_pdf, title="Saturation analysis")

for (i in 1:length(prefixes)) {
	prefix <- prefixes[i]
	full_file <- full_files[i]
	samp_name <- basename(prefix)
	
	pcts <- seq(from=10,to=90,by=10)
	thresholds <- c(0, 0.1, 0.5, 1)
	counts <- as.data.frame(matrix(0, nrow=length(thresholds), ncol=length(pcts)+1))
	colnames(counts) <- c(pcts, 100)
	rownames(counts) <- thresholds
	counts$threshold <- thresholds
	#counts <- data.frame(count=rep(0, length(pcts)*length(thresholds)), sample_pct=rep(pcts, length(thresholds)), threshold=rep(thresholds, each=length(pcts)))

	for (pct in pcts) {
		subset.data <- read.table(sprintf("%s.sub_%d.cufflinks_output/genes.fpkm_tracking", prefix, pct), header=T, sep="\t")
		subset.counts <- unlist(lapply(thresholds, function(x) length(which(subset.data$FPKM>x))))
		counts[as.character(thresholds), as.character(pct)] <- subset.counts
	}

	full.data <- read.table(full_file, header=T, sep="\t")
	full.counts <- unlist(lapply(thresholds, function(x) length(which(full.data$FPKM>x))))
	counts[as.character(thresholds), "100"] <- full.counts

	# melt
	df <- melt(counts, id.vars="threshold", variable.name="subset_pct")

	# plot
	g <- ggplot(data=df, aes(x=subset_pct, y=value, group=factor(threshold), color=factor(threshold))) + ylab("Number of expressed transcripts") + xlab("Read subset percentage") + geom_line() + geom_point() + ggtitle(sprintf("Saturation analysis (%s)", samp_name)) + ylim(c(min(df$value)*0.8, max(df$value)*1.2)) + scale_color_discrete(name ="FPKM threshold")
	print(g)
}

dev.off()





