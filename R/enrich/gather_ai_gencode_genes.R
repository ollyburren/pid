library(data.table)
library(GenomicRanges)

## read in data from Javierre et al.

g<-fread("/home/ob219/scratch/pid/Javierre_prioritised_genes.csv")
g.ai<-unique(subset(g,Trait %in% c('CD','CEL','MS','PBC','RA','SLE','T1D','UC'))$EnsemblID)
g.ht<-unique((subset(g,Trait=='HT')$EnsemblID))

## this code to exclude intersections between the two
both<-intersect(g.ai,g.ht)
g.ai<-g.ai[!g.ai %in% both] 
g.ht<-g.ht[!g.ht %in% both] 

## load in exonic regions
(load("/home/ob219/scratch/pid/gencode.v26.exons.gr.RData"))

## select AI gene exonic regions
ai.gr<-subset(e.gr,ensg %in% g.ai)
## there are some missing

## take the exonic regions forward as we want to extract these from the file
library(rtracklayer)
## here we could just filter to include those exons with a ccds ?
export(reduce(ai.gr),con="/home/ob219/scratch/pid/all_ai_exons.bed")
## peversely read back in to split by chromosome
o<-fread("/home/ob219/scratch/pid/all_ai_exons.bed")
lapply(split(o,o$V1),function(x){
	chr<-unique(x$V1)
	tmp<-x[,.(V1,V2,V3)]
	write.table(tmp,file=sprintf("/home/ob219/scratch/pid/all_ai_exons_by_chr/%s.bed",chr),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
})


ht.gr<-subset(e.gr,ensg %in% g.ht)
## take the exonic regions forward as we want to extract these from the file
library(rtracklayer)
## here we could just filter to include those exons with a ccds ?
export(reduce(ht.gr),con="/home/ob219/scratch/pid/all_ht_exons.bed")
## peversely read back in to split by chromosome
o<-fread("/home/ob219/scratch/pid/all_ht_exons.bed")
lapply(split(o,o$V1),function(x){
        chr<-unique(x$V1)
        tmp<-x[,.(V1,V2,V3)]
        write.table(tmp,file=sprintf("/home/ob219/scratch/pid/all_ht_exons_by_chr/%s.bed",chr),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
})

