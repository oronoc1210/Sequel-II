#!/usr/bin/env Rscript

args<-commandArgs(trailingOnly = T)
label <- args[1]
fslen <- args[2]
fclen <- args[3]


slen=read.table(fslen,sep='\t',header=F)
clen=read.table(fclen,sep='\t',header=F)

pdf(paste(label,'len_dist.pdf',sep='_'))
#plot(hist(slen$V1),xlab='subread length',main=label)
#ablines(density(slen$V1))
plot(density(slen$V1),xlab='subread length',main=label)
plot(density(clen$V1),xlab='CCS read length',main=label)
dev.off()

