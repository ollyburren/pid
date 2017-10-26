## code to process Ig GWAS summary statistics so we can run through CHICGP pipeline

library(data.table)
library(magrittr)

## read in Rafal's plink input

cvid.dir<-'/home/ob219/scratch/telomere/rafal_results'
fs<-list.files(path=cvid.dir,pattern='*.linear',full.names=TRUE)
all.res <- rbindlist(lapply(fs,fread))

# adjust for population inflation
all.res<-all.res[!is.na(all.res$P),]
all.res[,chi.sq:=STAT^2]
lambda.pid<-median(all.res$chi.sq)/0.456
all.res[,P.adj:=2 * pnorm(sqrt(chi.sq/lambda.pid),lower.tail=FALSE)]

## remove MHC
rm.idx<-with(all.res,which((between(BP,31e6,35e6) & CHR==6)))
all.res<-all.res[-rm.idx,]





#DATA.DIR<-'/home/ob219/scratch/Ig_gwas/raw'
PMI.OUT.DIR<-'/home/ob219/scratch/telomere/pmi'
GENE.SCORE.OUT.DIR<-'/home/ob219/scratch/telomere/gene_score'
SH.FILE<-'/home/ob219/scratch/cvid/computeGeneScores.sh'
#info.impute.thresh<-0.4
#maf.thresh<-0.01
## remove those with HWE p.val < chi2=25
#hwe.thresh<-2*pnorm(sqrt(25),lower.tail=FALSE)

#fs<-list.files(path=DATA.DIR,pattern="*.txt",full.names=TRUE)
## load in LD blocks
library(rtracklayer)
ld.gr<-import.bed("/home/ob219/scratch/DATA/JAVIERRE_GWAS/support/0.1cM_regions.b37.bed")



addLDBlock<-function(DT,ld.gr){
  dt.gr<-with(DT,GRanges(seqnames=Rle(chr),ranges=IRanges(start=pos,end=pos),idx=1:nrow(DT)))
  ol<-as.matrix(findOverlaps(dt.gr,ld.gr))
  DT[ol[,1],ld:=ol[,2]]
}

## robustly sum logs
logsum <- function(x) {
  my.max <- max(x) ##take out the maximum value in log form)
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}

Var.data.cc <- function(f, N, s) {
    1 / (2 * N * f * (1 - f) * s * (1 - s))
}

## compute approx bayes factors and resultant posterior probabilities
## based on the assumption of one causal variant in a region
approx.bf.p <- function(p,f, N, s,pi_i,type='CC') {
    if(type=="QUANT") {
      sd.prior <- 0.15
      V <- Var.data(f, N)
    } else {
      sd.prior <- 0.2
      V <- Var.data.cc(f, N, s)
    }
    z <- qnorm(0.5 * p, lower.tail = FALSE)
    ## Shrinkage factor: ratio of the prior variance to the total variance
    r <- sd.prior^2 / (sd.prior^2 + V)
    ## Approximate BF  # I want ln scale to compare in log natural scale with LR diff
    lABF = 0.5 * (log(1-r) + (r * z^2))
    ## tABF - to add one we create another element at the end of zero for which pi_i is 1
    tABF <- c(lABF,0)
    vpi_i<-c(rep(pi_i,length(lABF)),1)
    sBF <- logsum(tABF + log(vpi_i))
    exp(lABF+log(pi_i)-sBF)
}


## load in reference allele freq
(load("/home/ob219/scratch/as_basis/1KG_support/all_EUR_support.RData"))
all.eur[,pid:=paste(chr,position,sep=':')]
setkey(all.eur,pid)
all.eur[,maf:=a2.f]
all.eur[maf>0.5,maf:=1-maf]

processFile<-function(DT,N=5760){
  #DT<-subset(DT,hwe>hwe.thresh & maf>maf.thresh & info_impute > info.impute.thresh)
  setnames(DT,c('chr','rsid','pos','allele_a','test','NMISS','beta','stat','p','chi.sq','p.adj'))
  DT<-addLDBlock(DT,ld.gr)
  ## need to add MAF
  DT[,pid:=paste(chr,pos,sep=':')]
  setkey(DT,pid)
  DT.maf<-all.eur[DT]
  N<-ncases + ncontrols
  ## need to remove NA mafs
  DT.maf <- DT.maf[!is.na(maf),]
  DT.maf[,ppi:=approx.bf.p(p=p.adj,f=maf,N=N,pi_i=1e-4,type="QUANT"),by=ld]
  ## this can then be input into the main pipeline if we output to pmi format
  out.DT<-with(DT.maf,data.table(chr=chr,start=position-1,end=position,rsid=name,maf=maf,pval=p.adj,ppi=ppi,imp.snp.pos=0,imp.r2=1))
  options(scipen=999)
  pmi.o.file<-file.path(PMI.OUT.DIR,paste('telomere','pmi',sep='.'))
  write.table(out.DT,file=pmi.o.file,sep="\t",row.names=FALSE,quote=FALSE)
  options(scipen=0)
  ## attempt to create the correct shell command to run hierachical gene score
  int.data<-'/home/ob219/scratch/DATA/JAVIERRE_GWAS/RDATA/javierre_interactions.RData'
  frag.data<-'/home/ob219/scratch/DATA/JAVIERRE_GWAS/RDATA/javierre_frags.by.ld.RData'
  csnps.data<-'/home/ob219/scratch/DATA/JAVIERRE_GWAS/RDATA/javierre_csnps.by.ld.RData'
  sets.yaml<-'/home/ob219/scratch/DATA/JAVIERRE_GWAS/support/javierre_tree.yaml'
  cmd_string<-"Rscript /home/ob219/git//CHIGP//R/computeGeneScoreH.R --pmi_file=%s --out_dir=%s --int=%s --frags=%s --csnps=%s --target.gene.cSNPs.only=1 --sets=%s --include.interactions.only=1 --decompose=1 --ppi.thresh=0.01 --BF.thresh=3"
  sprintf(cmd_string,pmi.o.file,GENE.SCORE.OUT.DIR,int.data,frag.data,csnps.data,sets.yaml)
}

processFile(all.res)

## read in file and process
top.hits<-res[res$overall_ppi>0.5 & biotype=='protein_coding',]

library(biomaRt)

all.hits <- res[biotype=='protein_coding',]

ensembl_archive <- 'feb2014.archive.ensembl.org'
ensembl <- useMart(host=ensembl_archive, biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
gene.details<-getBM(c("ensembl_gene_id","gene_biotype","external_gene_id","hgnc_symbol", "chromosome_name","start_position", "end_position"),
  mart=ensembl,filters=list(ensembl_gene_id=all.hits$ensg))
gene.details<-data.table(gene.details)
setkey(gene.details,ensembl_gene_id)
setkey(all.hits,ensg)
all.hits<-all.hits[gene.details]
all.hits<-all.hits[order(chromosome_name,start_position,overall_ppi),]
all.hits<-all.hits[,c(1:3,7,11,16:18),with=FALSE]
setnames(all.hits,c('Ensembl_ID','Biotype','Gene_Name','Type','Score','Gene_Chr','Gene_Start','Gene_End'))
all.hits<-all.hits[order(Score,Gene_Chr,Gene_Start,decreasing=TRUE),]
write.table(all.hits,file="~/tmp/CVID.csv",sep=",",row.names=FALSE,quote=FALSE)
write.table(all.hits[Score>0.5,],file="~/tmp/CVID_gt_0.5.csv",sep=",",row.names=FALSE,quote=FALSE)


#MH plot

library(qqman)

## downsample the boring part
THRESHOLD<-1e-2

DT<-DT[order(DT$chr,DT$pos),]
tmp<-split(DT$pos,DT$chr)


cs<-c(0,head(cumsum(as.numeric(sapply(tmp,max))),-1))
for(i in seq_along(tmp)){
  tmp[[i]]<-tmp[[i]] + cs[i]
}

DT$gpos<-do.call('c',tmp)
interesting<-which(DT$p.adj<=THRESHOLD)
boring<-DT[p.adj>THRESHOLD,]
#keep.idx <- do.call('c',
#  lapply(split(1:nrow(boring),cut(boring$gpos,quantile(boring$gpos,seq(0,1,by=1e-4)))),sample,10))
toplot<-DT[interesting,.(rsid,chr,pos,p.adj)]
setnames(toplot,c('SNP','CHR','BP','P'))
manhattan(toplot)
dev.print(pdf,'~/tmp/CVID_mh.pdf')

## MY CODE TO DO THIS

if(FALSE){
library(ggplot2)

DT<-DT[order(DT$chr,DT$pos),]
tmp<-split(DT$pos,DT$chr)


cs<-c(0,head(cumsum(as.numeric(sapply(tmp,max))),-1))
for(i in seq_along(tmp)){
  tmp[[i]]<-tmp[[i]] + cs[i]
}

DT$gpos<-do.call('c',tmp)

## for manhatten plotting we don't care about p-values under a certain threshold
THRESHOLD<-1e-2
ggplot(DT[DT$p.adj<THRESHOLD],aes(x=gpos,y=-log10(p.adj),color=chr%%2==0)) + geom_point(size=0.1) + theme_bw() +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
        geom_hline(yintercept = c(-log10(5e-8),-log10(1e-5)),color=c('red','blue')) + coord_cartesian(ylim=c(0,25))

}
