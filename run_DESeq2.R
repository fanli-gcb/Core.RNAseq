#! /usr/bin/Rscript

library(methods)
library("DESeq2")
library("gplots")
library("RColorBrewer")
library("ggplot2")
library("pheatmap")

args <- commandArgs(T)
if (length(args) < 6) {
	cat("USAGE: ./run_DESeq2.R counts_file factors_file out_txt out_up_txt out_down_txt out_pdf\n")
	q()
}

counts_file <- args[1]
factors_file <- args[2]
out_txt <- args[3]
out_up_txt <- args[4]
out_down_txt <- args[5]
out_pdf <- args[6]

# read in data
counts <- read.table(counts_file, header=T, sep="\t", row.names=1)
rows_to_keep <- setdiff(rownames(counts), rownames(counts)[grep("^__", rownames(counts))])
counts <- counts[rows_to_keep,]
counts <- round(counts)
condition <- factor(scan(factors_file, what="c", sep=","))

########################################################################################
pdf(out_pdf, title="DESeq2 plots")

dds <- DESeqDataSetFromMatrix(countData=counts, colData=data.frame(condition, row.names=colnames(counts)), design= ~condition)
colData(dds)$condition <- factor(colData(dds)$condition, levels=levels(condition))
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]
resSig <- subset(res, padj<0.05)
resSig.up <- subset(resSig, log2FoldChange>0)
resSig.down <- subset(resSig, log2FoldChange<0)

write.table(as.data.frame(resSig), file=out_txt, quote=F, sep="\t")
write.table(rownames(resSig.up), file=out_up_txt, quote=F, sep="\t", row.names=F, col.names=F)
write.table(rownames(resSig.down), file=out_down_txt, quote=F, sep="\t", row.names=F, col.names=F)


# diagnostic plots
plotDispEsts(dds)
DESeq2::plotMA(res, main="DESeq2 MA plot")

vsd = varianceStabilizingTransformation( dds, blind=TRUE )

# sampled gene clustering
select = rownames(resSig)
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(assay(vsd)[select,], Colv=F, dendrogram="row", col = hmcol, trace="none", margin=c(10, 6))

df <- as.data.frame(colData(dds)[,c("condition")]); colnames(df) <- "condition"
# all sighits
select = rownames(resSig)
pheatmap(assay(vsd)[select,], cluster_rows=T, show_rownames=F, cluster_cols=FALSE, main="All sighits")
tmp <- assay(vsd)[select,]; tmp <- t(apply(tmp, 1, function(x) (x-mean(x))/sd(x)))
pheatmap(tmp, cluster_rows=T, show_rownames=F, cluster_cols=FALSE, main="all sighits (Z-transformed)")
# down
select = rownames(resSig.up)
pheatmap(assay(vsd)[select,], cluster_rows=T, show_rownames=F, cluster_cols=FALSE, main="Down sighits")
tmp <- assay(vsd)[select,]; tmp <- t(apply(tmp, 1, function(x) (x-mean(x))/sd(x)))
pheatmap(tmp, cluster_rows=T, show_rownames=F, cluster_cols=FALSE, main="Down sighits (Z-transformed)")
# up
select = rownames(resSig.down)
pheatmap(assay(vsd)[select,], cluster_rows=T, show_rownames=F, cluster_cols=FALSE, main="Up sighits")
tmp <- assay(vsd)[select,]; tmp <- t(apply(tmp, 1, function(x) (x-mean(x))/sd(x)))
pheatmap(tmp, cluster_rows=T, show_rownames=F, cluster_cols=FALSE, main="Up sighits (Z-transformed)")


# heatmap of sample-to-sample distances
dists = dist( t( assay(vsd) ) )
mat = as.matrix( dists )
#rownames(mat) = colnames(mat) = with(pData(ddsBlind), paste(condition, libType, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))

# correlation matrix
correlation <- cor(assay(vsd))
corrstr <- apply(correlation, c(1,2), function(x) sprintf("%.3g", x))
colors2 <- colorRampPalette(c("white","darkblue"))(150)
heatmap.2(correlation, col=colors2, scale="none", trace="none", main="correlation between samples", cellnote=corrstr, notecex=0.4, notecol="black")

# PCA
plotPCA(vsd, intgroup="condition") + geom_text(aes(label = name), size = 1.5, col="black", hjust=0, vjust=1)

dev.off()


