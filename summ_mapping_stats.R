#!/usr/bin/Rscript

# load required packages
library(ggplot2)
library(reshape2)
library(plyr)

args <- commandArgs(T)
if (length(args) < 3) {
	cat("USAGE: ./summ_mapping_stats.R dir out_txt_fn out_pdf_fn\n")
	q()
}

dir <- args[1]
out_txt_fn <- args[2]
out_pdf_fn <- args[3]

samples <- list.files(sprintf("%s", dir))

trim.whitespace <- function(x) {
	sub("\\s+$", "", sub("^\\s+", "", x))
}


data <- {}
for (sample_id in samples) {
	fn <- sprintf("%s/%s/align_summary.txt", dir, sample_id)
	if (!file.exists(fn)) {
		samples <- setdiff(samples, sample_id)
		next
	}
	tmp <- readLines(fn)
	
	total.left <- as.numeric(unlist(strsplit(trim.whitespace(tmp[2]), "\\s+"))[3])
	mapped.left <- as.numeric(unlist(strsplit(trim.whitespace(tmp[3]), "\\s+"))[3])
	total.right <- as.numeric(unlist(strsplit(trim.whitespace(tmp[6]), "\\s+"))[3])
	mapped.right <- as.numeric(unlist(strsplit(trim.whitespace(tmp[7]), "\\s+"))[3])
	mapped.pairs <- as.numeric(unlist(strsplit(trim.whitespace(tmp[11]), "\\s+"))[3])
	if (length(tmp)==14) {
		mapped.discordant <- as.numeric(unlist(strsplit(trim.whitespace(tmp[13]), "\\s+"))[1])
	} else {
		mapped.discordant <- 0
	}
	mapped.concordant <- mapped.pairs - mapped.discordant
	mapped.overall <- round(mean(c(mapped.left, mapped.right)))
	
	data <- rbind(data, c(mapped.left, mapped.right, mapped.overall, total.left, mapped.pairs-mapped.discordant))
}
colnames(data) <- c("left_mapped", "right_mapped", "overall_mapped", "total_reads", "concordant_mapped")
rownames(data) <- samples
data <- as.data.frame(data)
data$pct_concordant <- data$concordant_mapped / data$total_reads

pdf(out_pdf_fn, width=20)

# distribution of total read counts
breaks <- pretty(range(data$total_reads), n = nclass.FD(data$total_reads)*3, min.n = 1)
bwidth <- breaks[2]-breaks[1]
ggplot(data, aes(x=total_reads)) + geom_histogram(aes(y = ..density..), binwidth=bwidth) + geom_density() + ggtitle("Total read counts") + xlim(breaks[1], breaks[length(breaks)])
# distribution of left mapped
ggplot(data, aes(x=left_mapped)) + geom_histogram(aes(y = ..density..), binwidth=bwidth) + geom_density() + ggtitle("Left mapped counts") + xlim(breaks[1], breaks[length(breaks)])
# distribution of right mapped
ggplot(data, aes(x=right_mapped)) + geom_histogram(aes(y = ..density..), binwidth=bwidth) + geom_density() + ggtitle("Right mapped counts") + xlim(breaks[1], breaks[length(breaks)])
# distribution of overall mapped
ggplot(data, aes(x=overall_mapped)) + geom_histogram(aes(y = ..density..), binwidth=bwidth) + geom_density() + ggtitle("Overall mapped counts") + xlim(breaks[1], breaks[length(breaks)])
# percentage concordant
ggplot(data, aes(x=pct_concordant)) + geom_histogram(aes(y = ..density..)) + geom_density() + ggtitle("Pct concordant mapped")
# scatterplot - %concordant vs total read count
ggplot(data, aes(x=total_reads, y=pct_concordant)) + geom_point() + stat_smooth(method="lm")

# boxplots of all numbers
df <- melt(data[, c("left_mapped", "right_mapped", "overall_mapped", "total_reads", "concordant_mapped")])
ggplot(df, aes(x=variable, y=value)) + geom_boxplot() + theme_classic()

# individual barplots
data$ID <- rownames(data)
data2 <- data[, c("ID", "concordant_mapped")]
data2$overall_mapped <- data$overall_mapped - data$concordant_mapped
data2$unmapped <- data$total_reads - (data2$overall_mapped + data2$concordant_mapped)
data2$ID <- factor(data2$ID, levels=rownames(data)[order(data$total_reads, decreasing=T)])
df <- melt(data2)
colnames(df) <- c("ID", "readType", "numReads")
ggplot(df, aes(x=ID, y=numReads, fill=readType)) + geom_bar(stat="identity", aes(order=readType)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_brewer(palette="Set1")

dev.off()

write.table(data, file=out_txt_fn, quote=F, sep="\t", row.names=T, col.names=T)


