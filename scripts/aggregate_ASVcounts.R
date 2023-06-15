#!/usr/bin/Rscript

library(tidyverse)

args<-commandArgs(T)
# get the list of number/codes of each sequencing group

perc <- 0.005
# set percentage of low-abundance cut-off in respect to the total amount of reads in each sample

minabu <- 5
# set minimum abundance of each ASV

otu.tabs <- list()
# create new list that will have all the otu tables of each sequencing run defined by the "groups" file before the denoising

# The following for loop will read each denoised ASV table (produced by recount_swarm.R) and convert to 0 the abundance of
# those sequences with an abundance lower the "perc" in respect to the total reads of each sample.
# Each table is then stored in the same list

for (i in args) {
	
	count_tab <- read.table(paste("denoising/",i,"/denoise_aggr.counts.tsv", sep=""), header=T, check.names=T, sep="\t")
	# "check names" is set to true, as the first data item "#OTU ID" was removed in the main bash script
	
	samples <- ncol(count_tab)
	
	for (y in seq(from=2, to=samples)) {
		
		tot <- round((sum(count_tab[,y])*perc)/100)
		
		count_tab[,y][count_tab[,y] <= tot] <- 0
		
	}
	
	otu.tabs[[i]] <- count_tab
	
}

otu.tot <- as.data.frame(bind_rows(otu.tabs) %>% group_by(X) %>% summarise_all(., sum, na.rm = TRUE))
# collapse the entire otu table by ASV name and sum the abundances. The column with the sequences names is renamed X by the
# "read.table" command thanks to "check.names=T"

names(otu.tot)[1] <- ""
# remove the name of the columns with the sequence names, here "X"

rownames(otu.tot) <- otu.tot[,1]
# set the first column data items as rownames

samples <- ncol(otu.tot)
otu.tot <- otu.tot[,c(seq(from=2, to=samples))]
# remove the first column

otu.tot.rem <- otu.tot[!(rowSums(otu.tot) < minabu),]
# remove all rows with a rowSum lower then minabu

# print otu table
write.table(otu.tot.rem,"denoising/ASV_counts_total.tab",sep="\t",quote=FALSE,row.names=TRUE,col.names=NA)



