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
	DT<-fread(canvas[1])
	DT<-subset(DT,ILMN_ID %in% id.file$Illumina.ID)
	DT$id<-1:nrow(DT)
	with(DT,GRanges(seqnames=Rle(CANVAS_CHROM),ranges=IRanges(start=CANVAS_START,end=CANVAS_END),id=id))
}

all.canvas<-lapply(canvas,processCANVAS)
names(all.canvas)<-gsub('\\.txt','',basename(canvas))

load(file.path(data.dir,'hnisz_SE_annotation.RData'))
filterAndAddSEPCHiC<-function(m,sv.gr){
##load the pchic SE overlaps merge.results.se
	ol<-as.matrix(findOverlaps(m$se.gr,sv.gr))
	cbind(DT[ol[,2],],m[ol[,1],c('ensg','name','pchic','se.tissue')])
}
results.canvas<-lapply(all.canvas,function(sv.gr){
	filt.canvas<-lapply(merge.results.se,filterAndAddSEPCHiC,sv.gr=sv.gr)
})

## TODO 
## checking that the code is doing the correct thing.
## how to usefully sumarize.
## what is the background are SV enriched in PID cases vs the other cohorts.
## how much of SE is duplicated/deleted - what is the GT is it het or hom. what about qualities what about gene deletions as the mmore parsimonious explanation




manta<-mfiles[grepl("manta",mfiles)]
