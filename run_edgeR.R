#! /usr/bin/Rscript

library(methods)
library("edgeR")
library("gplots")
library("RColorBrewer")
library("ggplot2")
library("pheatmap")

args <- commandArgs(T)
if (length(args) < 6) {
	cat("USAGE: ./run_edgeR.R counts_file factors_file out_txt out_up_txt out_down_txt out_pdf\n")
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
group <- factor(scan(factors_file, what="c", sep=","))

########################################################################################

y <- DGEList(counts=counts, group=group)
y <- calcNormFactors(y)
keep <- rowSums(cpm(y)>1) >= 2 # expressed CPM>1 in at least two samples
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y, design, robust=T)

res <- exactTest(y)
res$table$padj <- p.adjust(res$table$PValue, method="fdr")
inds <- decideTestsDGE(res, adjust.method="BH", p.value=0.05)
resSig <- res[which(inds!=0), ]
resSig.up <- res[which(inds==1), ]
resSig.down <- res[which(inds==-1), ]

write.table(as.data.frame(resSig$table), file=out_txt, quote=F, sep="\t")
write.table(rownames(resSig.up$table), file=out_up_txt, quote=F, sep="\t", row.names=F, col.names=F)
write.table(rownames(resSig.down$table), file=out_down_txt, quote=F, sep="\t", row.names=F, col.names=F)

# diagnostic plots
pdf(out_pdf, title="edgeR plots")

colors <- c("red", "blue")
plotMDS(y, col=colors[group], pch=0, labels=colnames(y))
plotBCV(y)

plotSmear(y, de.tags=rownames(resSig))


dev.off()


