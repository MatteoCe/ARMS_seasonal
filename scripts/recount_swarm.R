#!/usr/bin/Rscript

# Script for recount abundances of a dataset clustered using Swarm to a tabulated tsv file
# The script will read two arguments from the command line: the input file name (output file from swarm, without the
# abundance annotation added to the zotus name) and the abundances count table file from --otutabout of
# vsearch --usearch_global.
# The output will be a ".tsv" file with recalculated abundances per MOTU and per sample.
# Modified from Owen S. Wangensteen - "owi_recount_swarm.R" in GitHub Project Metabarpark 2017

args<-commandArgs(T)

if (length(args)<2) {message("Please enter a cluster list output file obtained from SWARM and a tabulated counts file from obitab.")} else
{  
  fileswarm <- args[1] #Must be an output file from swarm 
  filetab <- args[2] 
  outfile <-paste(fileswarm,".counts.tsv",sep="")   
  
  num_char_id <- 14
   
  # Read cluster list database
  message("Reading swarm database...")
  swarm_db <- readLines(fileswarm)
  clusters <- strsplit(swarm_db," ")
  total_swarms <- length(clusters)

 id <- NULL
  for (i in 1:total_swarms) for (j in 1:length(clusters[[i]])) {
    clusters[[i]][[j]] <- substr(clusters[[i]][[j]],1,num_char_id)  
    id[i] <- clusters[[i]][1]
  }
  names(clusters) <- id

  # Read counts database and keep only the needed clusters
  message("Reading tabulated database. This could take a while...")
  db <- read.table(filetab,sep="\t",head=T)
  names(db)[1] <- "id"
  numseqs <- nrow(db)
  samples <- length(names(db)[substr(names(db),1,6)!="id"])
  db.total <- merge(data.frame(id),db,by="id")
  id <- db.total$id
  numclust <- nrow(db.total)
  
  for (fila in 1:numclust){
    head <- id[fila]
    tails <- unlist(clusters[names(clusters)==head])
    db.reduced <- db[db$id %in% tails,]
    suma <- colSums(db.reduced[,substr(names(db.total),1,6)!="id", drop=FALSE])
    # mcecchetto included "drop=FALSE" to include group of samples with only one sample
    db.total[fila,substr(names(db.total),1,6)!="id"] <- suma
   }
  names(db.total[substr(names(db.total),1,6)!="id"]) <- substr(names(db.total[substr(names(db.total),1,6)!="id"]),8,nchar(names(db.total[substr(names(db.total),1,6)!="id"])))
  names(db.total)[1] <- ""
  rownames(db.total) <- db.total[,1]
  db.total <- db.total[,2:ncol(db.total), drop=FALSE]
  # mcecchetto included "drop=FALSE" to include group of samples with only one sample
  write.table(db.total,outfile,sep="\t",quote=F,row.names=T, col.names=NA)
message("File ", outfile, " written")
}

