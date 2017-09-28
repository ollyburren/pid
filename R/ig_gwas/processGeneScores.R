## code to parse the gene score data for Ig GWAS
library(data.table)
GENE.SCORE.OUT.DIR<-'/home/ob219/scratch/Ig_gwas/gene_score'
fs<-list.files(path=GENE.SCORE.OUT.DIR,pattern="*prioritised.tab",full.names=TRUE)

processFilePrior<-function(f){
  DT<-fread(f)
  DT<-subset(DT,biotype=='protein_coding')
  DT<-DT[order(DT$overall_ppi,decreasing=TRUE),]
  DT[,.(ensg,name,overall_ppi,node,isLeaf,disease)]
}

all.prioritised.dat<-rbindlist(lapply(fs,processFilePrior))

subset(all.prioritised.dat,overall_ppi > 0.5)

## check cand genes from Jonsson et al.
##cgenes<-subset(all.dat,name %in% c('RUNX3','CYP26B1','FADS1','FADS2','ZMAT5','UQCR10','ANKRD55','IL6ST','SH2B2','STAT1','FCGR2B','KLF10','SIGLEC1','GSDMB','ZPBP2','ORMDL3','MIEN1','IKZF3','PTPN22','RNF168','TNFSF13','BCL2L11','TRAF3','IGHA1','IGHG2','PPTN7','HVCN1','KLRC1','KLRC2','KLRC3','IRF8','LGMN','REL','MUTYH','TESK2','TLR1'))
##cgenes<-cgenes[order(cgenes$overall_ppi,decreasing=TRUE),]

## get scores for all.

fs<-list.files(path=GENE.SCORE.OUT.DIR,pattern="*full.tab",full.names=TRUE)

processFileFull<-function(f){
  DT<-fread(f)
  DT<-subset(DT,biotype=='protein_coding')
  DT<-DT[order(DT$overall_gene_score,decreasing=TRUE),]
  DT[,.(ensg,baitChr,name,overall_gene_score,disease)]
}

all.full.dat<-lapply(fs,processFileFull)

##

require(xlsx)
chet<-read.xlsx("/home/ob219/scratch/pid/RESULTS/possible_CHET.xlsx",sheetName = "all")
chet.DT<-data.table(chet)
## add in hic scores

for(i in seq_along(all.full.dat)){
  x<-all.full.dat[[i]]
  tname<-unique(x$disease)
  #idx<-match(x$ensg,chet.DT$ensg)
  idx<-match(chet.DT$ensg,x$ensg)
  chet.DT[[tname]]<- 0
  chet.DT[[tname]]<-x[idx,]$overall_gene_score
  chet.DT[is.na(chet.DT[[tname]]),][[tname]]<-0
}

## compute a dummy overall score

chet.DT[,Ig_overall:=1-prod(c(1-LN_IgA_Conc,1-LN_IgG_total,1-LN_IgM_Conc)),by="NA."]
chet.DT<-chet.DT[order(chet.DT$Ig_overall,decreasing=TRUE),]

## next add in the pli scores

pli<-fread('/home/ob219/scratch/pid/support/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt')
pli.DT<-pli[,list(max.pli=signif(max(pLI),digits=2)),by=gene]
chet.DT$pli<-pli.DT[match(chet.DT$gene,pli.DT$gene),]$max.pli
## overall rank
#chet.DT[is.na(pli),pli:=0]
write.xlsx(chet.DT,file="/home/ob219/scratch/pid/RESULTS/possible_CHET_Ig_scores.xlsx",sheetName = "all")
