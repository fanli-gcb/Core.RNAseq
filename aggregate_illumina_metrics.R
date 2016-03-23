#!/usr/bin/Rscript

# load required packages
library(reshape2)
library(plyr)
library(ggplot2)

args <- commandArgs(T)
if (length(args) < 2) {
	cat("USAGE: ./aggregate_illumina_metrics.R [directories with Illumina run folders] out_txt out_pdf\n")
	q()
}

nargs <- length(args)
out_txt <- args[nargs-1]
out_pdf <- args[nargs]


df <- {}
for (i in 1:(nargs-2)) {
	path <- args[i]
	runs <- basename(list.dirs(path, recursive=F))
	for (run in runs) {
		tmp <- system(sprintf("python /Lab_Share/fanli/code/MiNGS/run_illuminate.py %s/%s/", path, run), intern=T)
		if (length(tmp) > 0) {
			# remove warnings from Illuminate output
			inds_to_remove <- grep("^\\[", tmp)
			if (length(inds_to_remove) > 0) {
				tmp <- tmp[setdiff(1:length(tmp), inds_to_remove)]
			}
			out <- colsplit(tmp, "\\t", c("variable", "value"))
			out$run <- run
			out$instrument <- basename(path)
			df <- rbind(df, out)
		}
	}
}

convert.magic <- function(obj,types) {
	for (i in 1:length(obj)){
		FUN <- switch(types[i],character = as.character, numeric = as.numeric, factor = as.factor)
		obj[,i] <- FUN(obj[,i])
	}
	obj
}

df2 <- dcast(df, run+instrument~variable)
types <- c(rep("character", 4), rep("numeric", 23), "character", "numeric", "numeric")
df2 <- convert.magic(df2, types)
df2$Mean_Cluster_Density_Kmm2 <- df2$Mean_Cluster_Density / 1000 # convert to K/mm2
df2 <- subset(df2, !is.na(Overall_Percent_Q30))

# plotting
pdf(out_fn, width=12)
p <- ggplot(df2, aes(x=Mean_Cluster_Density_Kmm2, y=Percent_PF_Clusters)) + geom_point() + facet_wrap(~instrument, scales="free_x") + ggtitle("Cluster Density vs. %PF") + xlab("Cluster Density (K/mm2)") + ylab("%Cluster PF")
print(p)
p <- ggplot(df2, aes(x=Mean_Cluster_Density_Kmm2, y=Overall_Percent_Q30)) + geom_point() + facet_wrap(~instrument, scales="free_x") + ggtitle("Cluster Density vs. %Q30") + xlab("Cluster Density (K/mm2)") + ylab("%Q30")
print(p)
p <- ggplot(df2, aes(x=Mean_Cluster_Density_Kmm2, y=Total_PF_Clusters)) + geom_point() + facet_wrap(~instrument, scales="free_y") + ggtitle("Cluster Density vs. Total PF Clusters") + xlab("Cluster Density (K/mm2)") + ylab("Total PF Clusters")
print(p)

# also individual detailed plots for each instrument
for (instr in unique(df2$instrument)) {
	df3 <- subset(df2, instrument==instr)
	xoffset <- max(range(df3$Mean_Cluster_Density_Kmm2)[2] - range(df3$Mean_Cluster_Density_Kmm2)[1], 1)*0.01
	p <- ggplot(df3, aes(x=Mean_Cluster_Density_Kmm2, y=Percent_PF_Clusters, group=Chemistry)) + geom_point() + geom_text(aes(label=Experiment_Name), size=2, nudge_x=xoffset, hjust=0) + ggtitle(sprintf("Cluster Density vs. %%PF (%s)", instr)) + xlab("Cluster Density (K/mm2)") + ylab("%Cluster PF")
	print(p)
	p <- ggplot(df3, aes(x=Mean_Cluster_Density_Kmm2, y=Overall_Percent_Q30, group=Chemistry)) + geom_point() + geom_text(aes(label=Experiment_Name), size=2, nudge_x=xoffset, hjust=0) + ggtitle(sprintf("Cluster Density vs. %%Q30 (%s)", instr)) + xlab("Cluster Density (K/mm2)") + ylab("%Q30")
	print(p)
	p <- ggplot(df3, aes(x=Mean_Cluster_Density_Kmm2, y=Total_PF_Clusters, group=Chemistry)) + geom_point() + geom_text(aes(label=Experiment_Name), size=2, nudge_x=xoffset, hjust=0) + ggtitle(sprintf("Cluster Density vs. Total PF Clusters (%s)", instr)) + xlab("Cluster Density (K/mm2)") + ylab("Total PF Clusters")
	print(p)
}

dev.off()

write.table(df2, file=out_txt, quote=F, sep="\t", row.names=T, col.names=F)


