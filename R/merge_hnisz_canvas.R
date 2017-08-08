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

convertSVToGR<-function(DT){
	with(DT,GRanges(seqnames=Rle(CANVAS_CHROM),ranges=IRanges(start=CANVAS_START,end=CANVAS_END),canvas.id=id))
}

all.canvas<-lapply(canvas,processCANVAS)
names(all.canvas)<-gsub('\\.txt','',basename(canvas))

load(file.path(data.dir,'hnisz_SE_annotation_all_merged.RData'))
## convert to genomic ranges object
se.pchic.gr<-with(merge.results.se,GRanges(seqnames=Rle(se.gr.seqnames),
	ranges=IRanges(start=se.gr.start,end=se.gr.end),uid=uid,prey.no=unique.prey,
	p.ensg=p.ensg,p.genenames=p.names,p.biotype=p.bt,p.known.pid=p.known.pid))

## examine loss
loss.gr<-convertSVToGR(all.canvas[['canvascalls_allsamples_loss_20161012.ann.filt']])
## can add in annotations like this
## how many samples have the same enhancer removed ?
loss.merge<-data.table(as.data.frame(mergeByOverlaps(loss.gr,se.pchic.gr)))
loss.merge$bridge_id<-all.canvas[['canvascalls_allsamples_loss_20161012.ann.filt']][loss.merge$canvas.id,c('BRIDGE_ID'),with=FALSE]
## which super enhancers are deleted in more than one individual ?
by.samp<-sapply(split(loss.merge$bridge_id,loss.merge$uid),length)
more.one.ind<-loss.merge[loss.merge$uid %in% names(by.samp[by.samp>1]),]
## of the above do not overlap known PID genes
more.one.ind[!grep('1',more.one.ind$p.known.pid),]
## link genes to individuals
ts<-split(loss.merge$p.ensg,loss.merge$bridge_id)
genes2ind<-data.table(do.call('rbind',lapply(seq_along(ts),function(x) cbind(names(ts)[x],do.call('c',strsplit(ts[[x]],','))))))
setnames(genes2ind,c('bridge_id','ensg'))
setkeyv(genes2ind,c('bridge_id','ensg'))
genes2ind[,list(ind.count=length(unique(bridge_id))),by=ensg]
genes2ind<-unique(genes2ind)
gene.ind.count<-genes2ind[,list(ind.count=length(unique(bridge_id))),by=ensg]
## take a gene centric view
## first create a list of unique genes with their details
all.genes<-data.table(do.call('cbind',lapply(c('p.ensg','p.genenames','p.biotype','p.known.pid'),function(c) do.call('c',strsplit(loss.merge[[c]],',')))))
setnames(all.genes,c('p.ensg','p.genenames','p.biotype','p.known.pid'))
setkey(all.genes,p.ensg)
all.genes<-unique(all.genes)
## merge in individual count data
setkey(gene.ind.count,ensg)
all.genes<-gene.ind.count[all.genes]
## sort by most genes with most deletions across the cohort
all.genes<-all.genes[order(all.genes$ind.count,decreasing=TRUE),]
write(all.genes[all.genes$ind.count>1,]$p.genenames,file=file.path(data.dir,'all_deleted_genes.txt'))

## examine gain

gain.gr<-convertSVToGR(all.canvas[['canvascalls_allsamples_gain_20161012.ann.filt']])
## can add in annotations like this
## how many samples have the same enhancer removed ?
gain.merge<-data.table(as.data.frame(mergeByOverlaps(gain.gr,se.pchic.gr)))
gain.merge$bridge_id<-all.canvas[['canvascalls_allsamples_gain_20161012.ann.filt']][gain.merge$canvas.id,c('BRIDGE_ID'),with=FALSE]
## which super enhancers are deleted in more than one individual ?
by.samp<-sapply(split(gain.merge$bridge_id,gain.merge$uid),length)
more.one.ind<-gain.merge[gain.merge$uid %in% names(by.samp[by.samp>1]),]
## of the above do not overlap known PID genes
more.one.ind[!grep('1',more.one.ind$p.known.pid),]
## link genes to individuals
ts<-split(gain.merge$p.ensg,gain.merge$bridge_id)
genes2ind<-data.table(do.call('rbind',lapply(seq_along(ts),function(x) cbind(names(ts)[x],do.call('c',strsplit(ts[[x]],','))))))
setnames(genes2ind,c('bridge_id','ensg'))
setkeyv(genes2ind,c('bridge_id','ensg'))
genes2ind[,list(ind.count=length(unique(bridge_id))),by=ensg]
genes2ind<-unique(genes2ind)
gene.ind.count<-genes2ind[,list(ind.count=length(unique(bridge_id))),by=ensg]
## take a gene centric view
## first create a list of unique genes with their details
all.genes<-data.table(do.call('cbind',lapply(c('p.ensg','p.genenames','p.biotype','p.known.pid'),function(c) do.call('c',strsplit(gain.merge[[c]],',')))))
setnames(all.genes,c('p.ensg','p.genenames','p.biotype','p.known.pid'))
setkey(all.genes,p.ensg)
all.genes<-unique(all.genes)
## merge in individual count data
setkey(gene.ind.count,ensg)
all.genes<-gene.ind.count[all.genes]
## sort by most genes with most deletions across the cohort
all.genes<-all.genes[order(all.genes$ind.count,decreasing=TRUE),]
write(all.genes[all.genes$p.biotype=='protein_coding',]$p.genenames,file=file.path(data.dir,'gain_genes.txt'))
