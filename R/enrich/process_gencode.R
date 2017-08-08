library(data.table)
library(GenomicRanges)

## gencode
# wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/GRCh37_mapping/gencode.v26lift37.annotation.gtf.gz

g<-fread("/home/ob219/scratch/pid/gencode.v26lift37.annotation.gtf")
g.exon<-subset(g,V3=='exon')
info<-do.call('rbind',lapply(strsplit(g.exon$V9,";"),function(f){
	gid<-sub('gene_id \"(.*)\"','\\1',f[grep("gene_id",f)])
	tid<-sub(' transcript_id \"(.*)\"','\\1',f[grep("transcript_id",f)])
	gt<-sub(' gene_type \"(.*)\"','\\1',f[grep("gene_type",f)])
	tt<-sub(' transcript_type \"(.*)\"','\\1',f[grep("transcript_type",f)])
	gn<-sub(' gene_name \"(.*)\"','\\1',f[grep("gene_name",f)])
	return(cbind(gid,tid,gt,tt,gn))
}))

final<-cbind(g.exon,info)

final.pc<-subset(final,gt=='protein_coding' & tt=='protein_coding')
final.pc$ensg<-sub("(ENSG[0-9]+)\\..*","\\1",final.pc$gid)
e.gr<-with(final.pc,GRanges(seqnames=Rle(sub("chr","",V1,fixed=TRUE)),ranges=IRanges(start=V4,end=V5),strand=V7,ensg=ensg,name=gn))
save(e.gr,file="/home/ob219/scratch/pid/gencode.v26.exons.gr.RData")

## how to compute tss ?
g.transcript<-subset(g,V3=='transcript')
info<-do.call('rbind',lapply(strsplit(g.transcript$V9,";"),function(f){
        gid<-sub('gene_id \"(.*)\"','\\1',f[grep("gene_id",f)])
        tid<-sub(' transcript_id \"(.*)\"','\\1',f[grep("transcript_id",f)])
        gt<-sub(' gene_type \"(.*)\"','\\1',f[grep("gene_type",f)])
        tt<-sub(' transcript_type \"(.*)\"','\\1',f[grep("transcript_type",f)])
        gn<-sub(' gene_name \"(.*)\"','\\1',f[grep("gene_name",f)])
        return(cbind(gid,tid,gt,tt,gn))
}))

final.t<-cbind(g.transcript,info)
final.tpc<-subset(final.t,gt=='protein_coding' & tt=='protein_coding')
final.tpc$ensg<-sub("(ENSG[0-9]+)\\..*","\\1",final.tpc$gid)
t.gr<-with(final.tpc,GRanges(seqnames=Rle(sub("chr","",V1,fixed=TRUE)),ranges=IRanges(start=V4,end=V5),strand=V7,ensg=ensg,name=gn))
save(t.gr,file="/home/ob219/scratch/pid/gencode.v26.transcript.gr.RData")

## get a list of tss

tss.gr<-t.gr
pos.ix<-which(strand(tss.gr)=='+')
neg.ix<-which(strand(tss.gr)=='-')
end(tss.gr[pos.ix,])<-start(tss.gr[pos.ix,])
start(tss.gr[neg.ix,])<-end(tss.gr[neg.ix,])
save(tss.gr,file="/home/ob219/scratch/pid/gencode.v26.tss.gr.RData")
