library(Gviz)
library(GenomicInteractions)
library(biomaRt)
library(data.table)
library(magrittr)

## plot a given region

## these are the parameters for the function.


# connect to biomart
e75.genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",  host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

## here is the merged Hnisz et al. data - currently missing CD20 data (needs to be added)
load('/Users/oliver/hpc_scratch/pid/hnisz_SE_annotation_all_merged.RData')
## make a handy object so we can compute overlaps
se.gr<-with(merge.results.se,GRanges(seqnames=Rle(paste0('chr',se.gr.seqnames)),ranges=IRanges(start=se.gr.start,end=se.gr.end),uid=uid))

tracks<-list()

## read in the interaction data and generate helper objects
hic<-fread("/Users/oliver/hpc_scratch/DATA/JAVIERRE_GWAS/chic/merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab")
baits.gr<-unique(hic[,.(baitChr,baitStart,baitEnd,baitID,ensg),key='baitID']) %>% with(.,GRanges(seqnames=Rle(paste0('chr',baitChr)),ranges=IRanges(start=baitStart,end=baitEnd),id=baitID,ensg=ensg))
pirs.gr<-unique(hic[,.(oeChr,oeStart,oeEnd,oeID,ensg),key='oeID']) %>% with(.,GRanges(seqnames=Rle(paste0('chr',oeChr)),ranges=IRanges(start=oeStart,end=oeEnd),id=oeID))

## what interaction do we want, it's specific in that it overlaps the SV and has a bait that is the same as the target gene

getInteraction<-function(gr){
    oeid<-subsetByOverlaps(pirs.gr,gr)$id
    bids<-baits.gr[baits.gr$ensg==gr$ensg,]$id
    fhic<-subset(hic,(oeID %in% oeid & baitID %in% bids) & biotype=='protein_coding')
    b.gr<-with(fhic,GRanges(seqnames=Rle(paste0('chr',baitChr)), IRanges(start=baitStart, end=baitEnd)))
    o.gr<-with(fhic,GRanges(seqnames=Rle(paste0('chr',oeChr)), IRanges(start=oeStart, end=oeEnd)))
    GenomicInteractions(b.gr, o.gr, counts=1)
    
}



getGenes<-function(interaction){
    intlength<-calculateDistances(interaction,method='outer')
    intStart<-start(anchorOne(interaction))
    intChr<-seqnames(anchorOne(interaction))
    genesToGet<-unique(subsetByOverlaps(baits.gr,GRanges(seqnames=intChr,ranges=IRanges(start=intStart,width=intlength)))$ensg)
    BiomartGeneRegionTrack(
        genome="hg19", name="Genes", transcriptAnnotation="symbol",
        mart=e75.genemart,
        collapseTranscripts="meta",shape = "arrow",
        stackHeight=0.2, filters=list(with_ox_refseq_mrna=T,ensembl_gene_id=genesToGet))
}

getSE<-function(interaction){
    pir.gr<-anchorTwo(interaction)
    subsetByOverlaps(se.gr,pir.gr)
}


plotter<-function(s.gr){
    tracks<-list()
    interaction<-getInteraction(s.gr)
    genes<-getGenes(interaction)
    se<-getSE(interaction)
    if(nrow(as.matrix(findOverlaps(s.gr,se))) ==0){
        sv<-AnnotationTrack(s.gr,name='SV',id=s.gr$individual,featureAnnotation="id",fontcolor.feature = "black",cex.feature=0.7)
    }else{
        sv<-AnnotationTrack(s.gr,name='SV',id=s.gr$individual,featureAnnotation="id",fontcolor.feature = "black",cex.feature=0.7,fill='red',col='red')
    }
    ## next add the tracks so we can plot
    ## axis first
    tracks$axis<-GenomeAxisTrack()
    ## interactions
    tracks$int<-InteractionTrack(name='Int', interaction,chromosome=as.character(seqnames(s.gr)))
    displayPars(tracks$int)$anchor.height <- 0.2
    displayPars(tracks$int)$col.anchors.line <- "darkgrey"
    displayPars(tracks$int)$col.anchors.fill <- "lightgrey"
    displayPars(tracks$int)$col.outside <- "white"
    displayPars(tracks$int)$col.interactions <- 'black'
    ## genes next
    tracks$gene<-genes
    ## structural variant
    tracks$sv<-sv
    ## super enhancers
    ## get the sample types
    seLabels<-sub("([^:]+):.*","\\1",se$uid)
    tracks$se<-AnnotationTrack(name='SE',se,id=seLabels, fontcolor.feature = "darkblue",
                               featureAnnotation="id", cex.feature = 0.7)
    paste(s.gr$ensg,s.gr$individual,sep=' ')
    plotTracks(tracks,main=paste(s.gr$ensg,s.gr$individual,sep=' '))
}

## can we plot for all ?

all.pse.chets<-fread('/Users/oliver/hpc_scratch/pid/RESULTS/possible_CHET.csv')
all.pse.chets<-subset(all.pse.chets,!is.na(sv.chr))
setkey(all.pse.chets,individual)
all.pse.chets<-unique(all.pse.chets)
## ok create the correct things 
pse.gr<-with(all.pse.chets,GRanges(seqnames=Rle(paste0('chr',sv.chr)),ranges=IRanges(start=sv.start,end=sv.end),ensg=target.gene,individual=individual))


pdf("~/tmp/priorGenes.pdf")
for(i in 1:length(pse.gr)){
    plotter(pse.gr[i,])   
}
dev.off()
