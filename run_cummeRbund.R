#!/usr/bin/Rscript

# load required packages
library(cummeRbund)

args <- commandArgs(T)
if (length(args) < 5) {
	cat("USAGE: ./run_cummeRbund.R cuffdiff_dir conditions transcripts_gtf_file genome_version out_pdf\n")
	q()
}

cuffdiff_dir <- args[1]
conditions_str <- args[2]
conditions <- unlist(strsplit(conditions_str, ","))
transcripts_gtf_file <- args[3]
genome_version <- args[4]
out_pdf <- args[5]

# Initialize cummeRbund database
cuff <- readCufflinks(dir = cuffdiff_dir, gtfFile = transcripts_gtf_file, genome = genome_version)
pdf(out_pdf, title=sprintf("cummeRbund output - %s", conditions_str))

# dispersion plot
disp <- dispersionPlot(genes(cuff))
disp

# CV2 plots
genes.scv <- fpkmSCVPlot(genes(cuff))
genes.scv
#isoforms.scv <- fpkmSCVPlot(isoforms(cuff))
#isoforms.scv

# FPKM density and boxplots
dens.gene <- csDensity(genes(cuff))
dens.gene
dens.gene.rep <- csDensity(genes(cuff),replicates=T)
dens.gene.rep
#dens.isoform <- csDensity(isoforms(cuff))
#dens.isoform
#dens.isoform.rep <- csDensity(isoforms(cuff),replicates=T)
#dens.isoform.rep
box.gene.rep <- csBoxplot(genes(cuff),replicates=T)
box.gene.rep
#box.isoform.rep <- csBoxplot(isoforms(cuff),replicates=T)
#box.isoform.rep

# pairwise log10(FPKM) scatterplots
s <- csScatterMatrix(genes(cuff))
s

# dendrograms
dend <- csDendro(genes(cuff))
dend
dend.rep <- csDendro(genes(cuff),replicates=T)
dend.rep

# MA plots
cond1 <- conditions[1]
cond2 <- conditions[2]
m <- MAplot(genes(cuff), cond1, cond2)

# volcano plots
v <- csVolcanoMatrix(genes(cuff))

# distance matrix
myDistHeat <- csDistHeat(genes(cuff))
myDistHeat
myDistHeat.rep <- csDistHeat(genes(cuff),replicates=T)
myDistHeat.rep

# PCA and MDS
genes.PCA <- PCAplot(genes(cuff),"PC1","PC2")
genes.PCA
genes.PCA.rep <- PCAplot(genes(cuff),"PC1","PC2",replicates=T)
genes.PCA.rep
#genes.MDS <- MDSplot(genes(cuff))
#genes.MDS
#genes.MDS.rep <- MDSplot(genes(cuff),replicates=T)
#genes.MDS.rep



# terminate device driver
dev.off()
