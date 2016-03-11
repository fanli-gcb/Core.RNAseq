#!/usr/bin/Rscript

# load required packages
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(grid)

args <- commandArgs(T)
if (length(args) < 2) {
	cat("USAGE: ./plot_RnaSeqMetrics.R metrics_file out_pdf\n")
	q()
}

fn <- args[1]
out_pdf <- args[2]

metrics <- read.table(fn, header=T, sep="\t", skip=5, nrows=1)
df <- melt(metrics[,11:15], variable.name="class")

labels <- unlist(lapply(1:length(df$class), function(i) sprintf("%s\n(%.2f%%)", gsub("PCT_","", as.character(df$class[i])), df$value[i]*100)))
y.breaks <- cumsum(df$value) - df$value/2

# plot
pdf(out_pdf, title="RNA-seq Metrics")

ggplot(data=df, aes(x=factor(1), y=value, fill=factor(class))) + geom_bar(width=1, stat="identity") + guides(fill=guide_legend(override.aes=list(color=NA))) + coord_polar(theta = "y") + theme(axis.text.y = element_blank(), axis.text.x=element_text(color='black'), axis.ticks = element_blank(), panel.grid  = element_blank(), axis.title = element_blank(), legend.position="none") + scale_y_continuous(breaks=y.breaks, labels=labels)

dev.off()





