library(data.table)
library(GenomicRanges)
library(xtable)
data.dir<-'/scratch/ob219/pid/'
## load in the file of bridge pid cases
id.file<-fread(file.path(data.dir,'BRIDGE-PID_20161012release_846_index_cases.csv'))
setnames(id.file,make.names(names(id.file)))

## we use Olga's latest CNV calls 
cnv<-fread("/scratch/WGS10K/data/WGS10K/us/SV-processing-V2all-upgrades/output/output_files/20170104/for_release/dels_strict/deletions.strict.allsamples.dedup.ann.filt.0.03.txt")

## filter to just contain the PID patients

pid.idx<-which(cnv$ILMN_ID %in% id.file$Illumina.ID)
pid.cnv<-cnv[pid.idx,]
other.cnv<-cnv[-pid.idx,]

## add a unique ID to pid.cnv
pid.cnv$uid<-1:nrow(pid.cnv)

## for more obvious ones create a list of hom cnv's 
canvas.gt.idx<-split(1:nrow(pid.cnv),pid.cnv$CANVAS_GT)
manta.gt.idx<-split(1:nrow(pid.cnv),pid.cnv$MANTA_GT)

## create genomic ranges objects based on manta and snv calls 

pid.manta.gr<-with(pid.cnv,GRanges(seqnames=Rle(MANTA_CHROM),ranges=IRanges(start=MANTA_START,end=MANTA_END),id=uid))
pid.canvas.gr<-with(pid.cnv,GRanges(seqnames=Rle(MANTA_CHROM),ranges=IRanges(start=CANVAS_START,end=CANVAS_END),id=uid))

load(file.path(data.dir,'hnisz_SE_annotation_all_merged.RData'))
## convert to genomic ranges object
se.pchic.gr<-with(merge.results.se,GRanges(seqnames=Rle(se.gr.seqnames),
	ranges=IRanges(start=se.gr.start,end=se.gr.end),uid=uid,prey.no=unique.prey,
	p.ensg=p.ensg,p.genenames=p.names,p.biotype=p.bt,p.known.pid=p.known.pid))

## how many samples have the same enhancer removed ?
manta.loss.merge<-data.table(as.data.frame(mergeByOverlaps(pid.manta.gr,se.pchic.gr)))
canvas.loss.merge<-data.table(as.data.frame(mergeByOverlaps(pid.canvas.gr,se.pchic.gr)))
## are any of these homozygous deletions if so which ?
table(manta.loss.merge$id %in% manta.gt.idx[['1/1']])
table(canvas.loss.merge$id %in% canvas.gt.idx[['1/1']])
## answer is no none of PID deletions that overlap super enhancer regions are homs
## which individuals have a deletion that takes out SE of a gene
manta.uid<-manta.loss.merge[grep(1,manta.loss.merge$p.known.pid),]$pid.manta.gr.id
manta.cnv.pid.gene<-subset(pid.cnv,uid %in% manta.uid)
manta.loss.no.gene<-manta.loss.merge[manta.loss.merge$pid.manta.gr.id %in% manta.cnv.pid.gene[is.na(manta.cnv.pid.gene$ensembl_gene_id_coll),]$uid,]
manta.loss.has.gene<-manta.loss.merge[manta.loss.merge$pid.manta.gr.id %in% manta.cnv.pid.gene[!is.na(manta.cnv.pid.gene$ensembl_gene_id_coll),]$uid,]
## for these can we remove tissue specificity and collapse by gene and individual to look further. For these individuals can we identify SNV that is impacts other preferably same PID gene in CHET type arrangement ?
## add individual

getGIndTab<-function(DT,cnv){
	## try and add in the se region so we can look it up backwards to work out distance from gene and tissues etc.
	tmp<-split(DT[,c('se.pchic.gr.seqnames','se.pchic.gr.start','se.pchic.gr.end','p.genenames','p.known.pid','p.ensg','p.biotype'),with=FALSE],cnv[DT[['pid.manta.gr.id']],][['BRIDGE_ID']])
	#tmp<-split(DT[,c('p.genenames','p.known.pid','p.ensg','p.biotype'),with=FALSE],cnv[DT[['pid.manta.gr.id']],][['BRIDGE_ID']])
	print(names(tmp))
	ret<-lapply(names(tmp),function(n){
		foo<-tmp[[n]]
		ret.dt<-lapply(names(foo),function(na){
			print(na)
			do.call('c',strsplit(as.character(foo[[na]]),","))
		})
		ret.dt<-data.table(do.call('cbind',ret.dt))
		setnames(ret.dt,names(foo))
		ret.dt$sample<-n
		ret.dt
	})
	ret<-rbindlist(ret)
	ret$k<-paste(ret$p.ensg,ret$sample,sep=':')
	setkey(ret,k)
	unique(ret)
}

unique.no.genes<-getGIndTab(manta.loss.no.gene,pid.cnv)
## next we should look these genes up to see if they have variants that disrupt the same genes 
## we should check first that the individual does not have another deletion taking out another PID gene.
unique.has.genes<-getGIndTab(manta.loss.has.gene,pid.cnv)
## add in a label as to whether gene is associated by just chromatin contacts or because of gene deletion
all.genes<-rbind(unique.no.genes,unique.has.genes)
all.genes$SE.only<-FALSE
all.genes[1:nrow(unique.no.genes),]$SE.only<-TRUE
## doesn't look as if there are two hits on the same gene one taking out coding and the other taking out enhancer
sapply(split(all.genes$SE.only,all.genes$k),length)
## however to do this properly will need to get the annotations for which exons overlap rather than relying on the SE data above 
## as this could be a long way removed from the gene !

## TODO get the genes annotated in the CNV file and repeat the above analysis

## do any of the subjects have insults to two pid genes ?
sapply(split(all.genes$p.ensg,all.genes$sample),length)
## some have 7 genes ! Probable artifact
## need to check these individuals to see if there is evidence for rare/common SNV in these genes that mean that gene is knocked out.

#library(jsonlite)

## for the person with 7 genes are these spatially close together

#getBaitCoords<-function(ensg){
#	url.gene<-sprintf("https://www.chicp.org/chicp/search?searchTerm=%s&targetIdx=CHICAGO",ensg)
#	message(url.gene)
#	chic.data.gene<-fromJSON(url.gene)
#	tmp<-unique(chic.data.gene$hic[,c('baitChr','baitStart_ori','baitEnd_ori')])
#	tmp$ensg<-ensg
#	tmp
#}

#baits.coord<-lapply(all.pid.genes[all.pid.genes$SE.only==1,]$p.ensg,getBaitCoords)
#baits.coord<-rbindlist(baits.coord)

## export as a latex table
#library(xtable)

#setkey(baits.coord,ensg)
#setkey(all.pid.genes,p.ensg)

#out.tab<-all.pid.genes[baits.coord]
out.tab<-all.genes
setkey(out.tab,k)
out.tab<-unique(out.tab)
#out.tab<-out.tab[order(baitChr,baitStart_ori,sample),]
out.tab<-out.tab[order(se.pchic.gr.seqnames,se.pchic.gr.start,sample),]
xtable(out.tab[,c('p.genenames','sample','SE.only','se.pchic.gr.seqnames','se.pchic.gr.start','se.pchic.gr.end'),with=FALSE])
## TOMORROW NEED TO COMPUTE THE OVERLAP BETWEEN se region and the deletion and add this to the file.
save(out.tab,file="/scratch/ob219/pid/loss_merge_genes/all.se.0.03.RData")
stop()
## do for others
## this code gets where the bait is - it's a bit slow and I think we can do better
baits.coord<-lapply(all.pid.genes[all.pid.genes$SE.only==0,]$p.ensg,getBaitCoords)
baits.coord<-rbindlist(baits.coord)

## export as a latex table
library(xtable)

setkey(baits.coord,ensg)
setkey(all.pid.genes,p.ensg)

out.tab<-all.pid.genes[baits.coord]
setkey(out.tab,k)
out.tab<-unique(out.tab)
out.tab<-out.tab[order(baitChr,baitStart_ori,sample),]
xtable(out.tab[,c('p.genenames','sample','SE.only','baitChr','baitStart_ori','baitEnd_ori'),with=FALSE])


stop()


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
