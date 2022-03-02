#!/usr/bin/env Rscript

args <- commandArgs(T);

ContigSize <- args[1]
MappedReadCount <- args[2]
RPKMoutput <- args[3]

ContigLength<-read.table(ContigSize, header=F)
MappedReads<-read.table(MappedReadCount, header=F)
Counts<-read.table(AbsoluteReadsCount, header=F)
rownames(Counts) <- as.character(Counts[, 1])
Counts <- Counts[, -1, drop = F]
ContigLength2<-do.call("cbind", replicate(nrow(MappedReads),ContigLength,simplify = FALSE))
MappedReads2<-do.call("rbind", replicate(nrow(ContigLength),as.data.frame(t(MappedReads)),simplify = FALSE))
MapXLength<-((ContigLength2/1000)*(MappedReads2/1000000))
RPKM<-(Counts/MapXLength)
write.table(RPKM,file=RPKMoutput,sep="\t", col.names = F, row.names = T)

