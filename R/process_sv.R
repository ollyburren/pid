library(data.table)
library(GenomicRanges)
pid.sv.dir<-'/scratch/WGS10K/data/release/20161012-A/merged-sv/all'
data.dir<-'/scratch/ob219/pid/'
mfiles<-list.files(path=pid.sv.dir,pattern='filt*.*txt',full.names=TRUE)
## load in the file of bridge pid cases
id.file<-fread(file.path(data.dir,'BRIDGE-PID_20161012release_846_index_cases.csv'))
setnames(id.file,make.names(names(id.file)))
##CANVAS - insertion deletions based on read depth
##MANTA - translocations, deletions, tandem duplications, insertions, and inversions based on both paired read fragment spanning and split read evidence
canvas<-mfiles[grepl("canvas",mfiles)]
## process canvas
processCANVAS<-function(f){
	DT<-fread(f)
	DT<-subset(DT,ILMN_ID %in% id.file$Illumina.ID)
	DT$id<-1:nrow(DT)
	DT
}

all.canvas<-lapply(canvas,processCANVAS)
names(all.canvas)<-gsub('\\.txt','',basename(canvas))

load(file.path(data.dir,'hnisz_SE_annotation.RData'))
filterAndAddSEPCHiC<-function(m,sv.gr,sv.dt){
##load the pchic SE overlaps merge.results.se
	ol<-as.matrix(findOverlaps(m$se.gr,sv.gr))
	cbind(sv.dt[ol[,2],],m[ol[,1],c('ensg','name','pchic','se.tissue')])
}

results.canvas<-lapply(all.canvas,function(sv){
	sv.gr<-with(sv,GRanges(seqnames=Rle(CANVAS_CHROM),ranges=IRanges(start=CANVAS_START,end=CANVAS_END),id=id))
	filt.canvas<-lapply(merge.results.se,function(mse)filterAndAddSEPCHiC(mse,sv.gr,sv))
})

## TODO 
## checking that the code is doing the correct thing.
## how to usefully sumarize.
## what is the background are SV enriched in PID cases vs the other cohorts.
## how much of SE is duplicated/deleted - what is the GT is it het or hom. what about qualities what about gene deletions as the mmore parsimonious explanation

## prototype downstream analysis on non exonic loss.
ft<-'canvascalls_allsamples_loss_20161012.ann.filt'
all.se.ft<-results.canvas[[ft]]
se.t.types<-sapply(merge.results.se,function(n) sub("\\.csv","",basename(unique(n$se.tissue))))
rows<-cbind(se.t.types,sapply(all.se.ft,nrow))
peeps<-cbind(se.t.types,sapply(all.se.ft,function(x) paste(unique(x$BRIDGE_ID),collapse=',')))
genes<-cbind(se.t.types,sapply(all.se.ft,function(x) paste(unique(sprintf("%s(%s)",x$name,x$BRIDGE_ID)),collapse=',')))
## all genes
sapply(all.se.ft,'[[','name')

manta<-mfiles[grepl("manta",mfiles)]
