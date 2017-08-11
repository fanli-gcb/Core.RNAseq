#! /usr/bin/Rscript

library(methods)
library("DESeq2")
library("gplots")
library("RColorBrewer")
library("reshape2")

plotPCAWithSampleNames = function (x, intgroup = "condition", ntop = 500, returnData = FALSE, select = NULL) 
{
	library("RColorBrewer")
	library("genefilter")
	library("ggplot2")
	rv <- rowVars(assay(x))
	if (is.null(select)) {
		select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
	}
	pca <- prcomp(t(assay(x)[select, ]))
	percentVar <- pca$sdev^2/sum(pca$sdev^2)
	if (!all(intgroup %in% names(colData(x)))) {
	  stop("the argument 'intgroup' should specify columns of colData(dds)")
	}
	intgroup.df <- as.data.frame(colData(x)[, intgroup, drop = FALSE])
	group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
	d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, intgroup.df, names = colnames(x))
	if (returnData) {
	  attr(d, "percentVar") <- percentVar[1:2]
	  return(d)
	}
	ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group", label = "names")) + geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + geom_text(size=3, hjust=0, vjust=0)
}

getPCALoadings = function (x,  ntop = 500, plottop=20, PC=1, select = NULL) 
{
	library("RColorBrewer")
	library("genefilter")
	library("ggplot2")
	rv <- rowVars(assay(x))
	if (is.null(select)) {
		select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
	}
	pca <- prcomp(t(assay(x)[select, ]))
	d <- pca$rotation[,PC]; d <- d[order(abs(d), decreasing=T)]
	d <- d[1:plottop]
	return(d)
}

mapping_fn <- "/Lab_Share/${PROJECT}/SAMPLE_SHEET.txt"
counts_dir <- "/Lab_Share/${PROJECT}/featureCounts"
out_pdf <- "/Lab_Share/${PROJECT}/DESeq2/DESeq2_results.pdf"

# read in data - featureCounts
mapping <- read.table(mapping_fn, header=T, as.is=T, sep="\t", row.names=1)
counts <- {}
for (sid in rownames(mapping)) {
	fn <- sprintf("%s/%s.featureCounts.counts.txt", counts_dir, sid)
	tmp <- read.table(fn, header=T, as.is=T, sep="\t", row.names=1)
	counts <- cbind(counts, tmp[,6])
}
gene_names <- rownames(tmp)
rownames(counts) <- gene_names
colnames(counts) <- rownames(mapping)
counts <- as.data.frame(counts)

# map to gene names
gene2name <- read.table("/Lab_Share/Ensembl/GRCh38.p7-85/genes2names.txt", header=F, as.is=T, sep="\t", row.names=1)
counts$gene_name <- gene2name[rownames(counts), "V2"]
counts <- aggregate(.~gene_name, counts, sum)
rownames(counts) <- counts$gene_name
counts <- counts[,-1]

# filter out genes w/ < 10 reads total
inds <- rowSums(counts) > 10
counts <- counts[inds,]

set.seed(sum(dim(counts)))

pdf(out_pdf, title="DESeq2 plots")

###########################################################################################
## overall PCA and diagnostic plots
dds <- DESeqDataSetFromMatrix(countData=counts, colData=mapping, design= ~Condition + CellLine)
dds <- DESeq(dds, parallel=T)

# k-means clustering
vsd = varianceStabilizingTransformation( dds, blind=TRUE )
out <- as.data.frame(assay(vsd))

# heatmap of sample-to-sample distances
dists = dist( t( assay(vsd) ) )
mat = as.matrix( dists )
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
#rownames(mat) = colnames(mat) = with(pData(ddsBlind), paste(condition, libType, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))

# correlation matrix
correlation <- cor(assay(vsd))
corrstr <- apply(correlation, c(1,2), function(x) sprintf("%.3g", x))
colors2 <- colorRampPalette(c("white","darkblue"))(150)
heatmap.2(correlation, col=colors2, scale="none", trace="none", main="correlation between samples", cellnote=corrstr, notecex=0.4, notecol="black")

# PCA
plotPCAWithSampleNames(vsd, intgroup="Condition")
plotPCA(vsd, intgroup="Condition")
plotPCAWithSampleNames(vsd, intgroup="CellLine")
plotPCA(vsd, intgroup="CellLine")
for (princ in c(1,2)) {
	d <- getPCALoadings(vsd, PC=princ)
	d <- data.frame(gene=names(d), loading=d)
	d <- d[order(d$loading),]
	d$gene <- factor(as.character(d$gene), levels=as.character(d$gene))
	p <- ggplot(d, aes(x=gene, y=loading)) + geom_bar(stat="identity") + scale_fill_manual(values="blue") + coord_flip() + theme_classic() + ggtitle(sprintf("PC %d", princ))
	print(p)
}

# differential analysis: by Condition
res <- results(dds, contrast=c("Condition", "Disease", "Control"))
res <- res[order(res$padj),]
res$gene_symbol <- gene_to_name[rownames(res), "Name"]
resSig <- subset(res, padj<0.05)
write.table(res, file="/Lab_Share/${PROJECT}/DESeq2/DESeq_results.Condition.txt", quote=F, sep="\t", row.names=T, col.names=T)


dev.off()





