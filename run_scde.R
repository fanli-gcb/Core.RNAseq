#!/usr/bin/Rscript

# load required packages
library(scde)

args <- commandArgs(T)
if (length(args) < 4) {
	cat("USAGE: ./run_scde.R counts_file factors_file out_txt out_pdf\n")
	q()
}

counts_file <- args[1]
factors_file <- args[2]
out_txt <- args[3]
out_pdf <- args[4]

# read in 
counts <- read.table(counts_file, header=T, sep="\t", row.names=1)
sg <- factor(scan(factors_file, what="c", sep=","))
names(sg) <- colnames(counts)

pdf(out_pdf, title="SCDE analysis")

## DE analysis
# omit genes that are never detected
cd <- counts
cd <- cd[rowSums(cd)>0,];
# omit cells with very poor coverage
cd <- cd[,colSums(cd)>1e4]; 

# number of local process cores to use in processing
n.cores <- 2;
# calculate models
o.ifm <- scde.error.models(counts=cd,groups=sg,n.cores=n.cores,threshold.segmentation=T,save.crossfit.plots=T,save.model.plots=T,verbose=1);
valid.cells <- o.ifm$corr.a >0;
table(valid.cells)
o.ifm <- o.ifm[valid.cells,];
o.prior <- scde.expression.prior(models=o.ifm,counts=cd,length.out=400,show.plot=T)
# define two groups of cells
groups <- sg
names(groups) <- row.names(o.ifm);
# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm,cd,o.prior,groups=groups,n.randomizations=100,n.cores=n.cores,verbose=1)

write.table(ediff[order(abs(ediff$Z),decreasing=T),],file=out_txt,row.names=T,col.names=T,sep="\t",quote=F)

## distance measure
p.self.fail <- scde.failure.probability(models=o.ifm,counts=cd)
# simulate drop-outs
# note: using 10 sampling rounds for illustration here. ~500 or more should be used.
n.simulations <- 10; k <- 0.9;
cell.names <- colnames(cd); names(cell.names) <- cell.names;
dl <- mclapply(1:n.simulations,function(i) {
  scd1 <- do.call(cbind,lapply(cell.names,function(nam) {
    x <- cd[,nam];
    # replace predicted drop outs with NAs
    x[!as.logical(rbinom(length(x),1,1-p.self.fail[,nam]*k))] <- NA;
    x;
    }))
  rownames(scd1) <- rownames(cd); 
  # calculate correlation on the complete observation pairs
  cor(log10(scd1+1),use="pairwise.complete.obs");
},mc.cores=n.cores)
# calculate average distance across sampling rounds
direct.dist <- as.dist(1-Reduce("+",dl)/length(dl))

pr <- prcomp(direct.dist)
samps <- sg[rownames(as.matrix(direct.dist))]
plot(pr$rotation[,"PC1"], pr$rotation[,"PC2"], col=c("red","darkgreen")[samps], pch=20, xlab="PC1", ylab="PC2", main="PCA of cells")
legend("bottomright", legend=levels(sg), fill=c("red","darkgreen"))

dev.off()


