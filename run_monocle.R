#!/usr/bin/Rscript

# load required packages
library(monocle)
library(reshape2)
library(ggplot2)

args <- commandArgs(T)
if (length(args) < 4) {
	cat("USAGE: ./run_monocle.R fpkm_table_file sample_sheet_file gene_attr_file out_pdf\n")
	q()
}

fpkm_table_file <- args[1]
sample_sheet_file <- args[2]
gene_attr_file <- args[3]
out_pdf <- args[4]

# read in data and set up objects
fpkm_matrix <- read.table(fpkm_table_file, header=T, sep="\t", row.names=1)
sample_sheet <- read.table(sample_sheet_file, header=T, sep="\t", row.names=1)
gene_ann <- read.table(gene_attr_file, header=T, row.names=1, sep="\t")
pd <- new("AnnotatedDataFrame", data=sample_sheet)
fd <- new("AnnotatedDataFrame", data=gene_ann)
cds <- newCellDataSet(as.matrix(fpkm_matrix), phenoData = pd, featureData = fd)

pdf(out_pdf, title="Monocle analysis")

# standard differential analysis
min_expression <- 0.1
cds <- detectGenes(cds, min_expr=min_expression)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 2))

samples_by_group <- aggregate(rownames(pData(cds)), by=list(as.character(pData(cds)$group)), FUN=paste)$x
tmp <- lapply(samples_by_group, function(x) apply(exprs(cds[, unlist(x)]), 1, function(r) any(r >= min_expression)))
expressed_in_all_groups <- unlist(apply(data.frame(matrix(unlist(tmp), nrow=length(unlist(tmp[1])), byrow=F)), 1, function(x) all(x)))
expressed_in_all_groups.genes <- rownames(exprs(cds))[expressed_in_all_groups]

L <- log(exprs(cds))
melted_dens_df <- melt(t(scale(t(L))))
qplot(value, geom = "density", data = melted_dens_df) + stat_function(fun = dnorm, size = 0.5, color = "red") + xlab("Standardized log(FPKM)") + ylab("Density")

diff_test_res <- differentialGeneTest(cds[expressed_in_all_groups.genes,], fullModelFormulaStr = "expression~group")
sig_genes <- subset(diff_test_res, qval < 0.1)
sig_genes <- merge(fData(cds), sig_genes, by="row.names")
rownames(sig_genes) <- sig_genes$Row.names
sig_genes <- sig_genes[,-1]

# ordering analysis
ordering_genes <- rownames(sig_genes)
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, use_irlba=F)
cds <- orderCells(cds, num_paths=1, reverse=F)
plot_spanning_tree(cds)

# genes that distinguish cell state
diff_test_res.state <- differentialGeneTest(cds[expressed_in_all_groups.genes,], fullModelFormulaStr = "expression~State")
diff_test_res.state <- merge(fData(cds), diff_test_res.state, by="row.names")

# genes that change as a function of pseudotime
diff_test_res.pt <- differentialGeneTest(cds[expressed_in_all_groups.genes,], fullModelFormulaStr = "expression~bs(Pseudotime)")
diff_test_res.pt <- merge(fData(cds), diff_test_res.pt, by="row.names")

# (optional) pseudotime plots of significant genes
my_genes <- as.character(subset(diff_test_res.pt, use_for_ordering==TRUE)$gene_id)
cds.subset <- cds[my_genes,]
for (gene in my_genes) {
	plot_genes_in_pseudotime(cds.subset[gene,], color_by="group")
}
# (optional) multi-factorial differential analysis

# clustering by pseudotime
full_model_fits <- fitModel(cds[expressed_in_all_groups.genes,], modelFormulaStr = "expression~bs(Pseudotime)")
expression_curve_matrix <- responseMatrix(full_model_fits)
clusters <- clusterGenes(expression_curve_matrix, k=4) # cluster::pam (partitioning around medoids) clustering with k=4
plot_clusters(cds[ordering_genes,], clusters)

#	(optional) violin plots of genes
for (gene in my_genes) {
	df <- melt(exprs(cds)[gene,])
	p <- ggplot(df, aes(x=1, y=value, fill="red")) + geom_violin() + ggtitle(sprintf("%s",gene))
	print(p)
}

dev.off()
