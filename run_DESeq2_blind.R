#! /usr/bin/Rscript

library(methods)
library("DESeq2")
library("gplots")
library("RColorBrewer")

args <- commandArgs(T)
if (length(args) < 3) {
	cat("USAGE: ./run_DESeq2_blind.R counts_file sample_sheet_file out_pdf\n")
	q()
}

counts_file <- args[1]
sample_sheet_file <- args[2]
out_pdf <- args[3]

# read in data
counts <- read.table(counts_file, header=T, sep="\t", row.names=1, strip.white=T)
rows_to_keep <- setdiff(rownames(counts), rownames(counts)[grep("^__", rownames(counts))])
counts <- counts[rows_to_keep,]
counts <- round(counts)
sample_sheet <- read.table(sample_sheet_file, header=T, sep="\t", row.names=1)
condition <- factor(sample_sheet[colnames(counts), "c1_run"])

########################################################################################
pdf(out_pdf, title="DESeq2 (blind) plots")

dds <- DESeqDataSetFromMatrix(countData=counts, colData=data.frame(condition, row.names=colnames(counts)), design= ~condition)
colData(dds)$condition <- factor(colData(dds)$condition, levels=levels(condition))

vsd = varianceStabilizingTransformation( dds, blind=TRUE )

# heatmap of sample-to-sample distances
dists = dist( t( assay(vsd) ) )
mat = as.matrix( dists )
#rownames(mat) = colnames(mat) = with(pData(ddsBlind), paste(condition, libType, sep=" : "))
#heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))

# correlation matrix
correlation <- cor(assay(vsd))
corrstr <- apply(correlation, c(1,2), function(x) sprintf("%.3g", x))
colors2 <- colorRampPalette(c("white","darkblue"))(150)
heatmap.2(correlation, col=colors2, scale="none", trace="none", main="correlation between samples", cellnote=corrstr, notecex=0.4, notecol="black")

# PCA
plotPCA(vsd, intgroup="condition")

dev.off()


