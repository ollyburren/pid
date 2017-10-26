library(Gviz)
library(GenomicInteractions)
library(biomaRt)
library(data.table)
library(magrittr)

## plot a given region

## these are the parameters for the function.


# connect to biomart
e75.genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",  host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

## hnisz
load('/Users/oliver/hpc_scratch/pid/hnisz_SE_annotation_all_merged.RData')
## make a handy object so we can compute overlaps
se.gr<-with(merge.results.se,GRanges(seqnames=Rle(paste0('chr',se.gr.seqnames)),ranges=IRanges(start=se.gr.start,end=se.gr.end),uid=uid))

## load in super enhancer data from blue print
## make a handy object so we can compute overlaps
se.bp.gr<-readRDS('/Users/oliver/hpc_scratch/pid/bp_se.RDS')
seqlevels(se.bp.gr)<-paste0('chr',seqlevels(se.bp.gr))

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

getSE<-function(interaction,gr){
    pir.gr<-anchorTwo(interaction)
    subsetByOverlaps(gr,pir.gr)
}

library(jsonlite)
library(httr)
ensg2symbol<-function(ensg){
  server <- "https://rest.ensembl.org"
  ext <- sprintf("/xrefs/id/%s?",ensg)
  request<-GET(paste(server, ext, sep = ""), content_type("application/json"))
  stop_for_status(request)
  jd<-fromJSON(toJSON(content(request)))
  unlist(subset(jd,dbname=='HGNC')[['display_id']])
}


plotter<-function(s.gr){
    tracks<-list()
    interaction<-getInteraction(s.gr)
    genes<-getGenes(interaction)
    se.bp<-getSE(interaction,se.bp.gr)
    se.hnisz<-getSE(interaction,se.gr)
    ## all se
    if(length(se.bp)!=0 & length(se.hnisz)!=0){
      all.se<-c(se.bp,se.hnisz)
    }else if (length(se.bp)==0 & length(se.hnisz)!=0){
      all.se <- se.hnisz
    }else if (length(se.bp)!=0 & length(se.hnisz)==0){
      all.se <- se.bp
    }
    if(nrow(as.matrix(findOverlaps(s.gr,all.se))) ==0){
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
    seLabels<-sub("([^:]+):.*","\\1",se.hnisz$uid)
    tracks$se_hnisz<-AnnotationTrack(name='SE Hnisz',se.hnisz,id=seLabels, fontcolor.feature = "darkblue",featureAnnotation="id", cex.feature = 0.7)
    tracks$se_bp<-AnnotationTrack(name='SE BP',se.bp,id=se.bp$uid, fontcolor.feature = "black",fill='orange',featureAnnotation="id", cex.feature = 0.7)
    plotTracks(tracks,main=paste(ensg2symbol(s.gr$ensg),s.gr$individual,sep=' '))
}

## can we plot for all ?

all.pse.chets<-fread('/Users/oliver/hpc_scratch/pid/RESULTS/PIK3C2B.csv')
all.pse.chets<-subset(all.pse.chets,!is.na(sv.chr))
setkey(all.pse.chets,individual)
all.pse.chets<-unique(all.pse.chets)
## ok create the correct things
pse.gr<-with(all.pse.chets,GRanges(seqnames=Rle(paste0('chr',sv.chr)),ranges=IRanges(start=sv.start,end=sv.end),ensg=target.gene,individual=individual))

## little function to convert ensg to gene symbol


pdf("~/tmp/PIK3C2B.pdf")
#pdf("~/tmp/priorGenes.pdf")
for(i in 1:length(pse.gr)){
    plotter(pse.gr[i,])
}
dev.off()
