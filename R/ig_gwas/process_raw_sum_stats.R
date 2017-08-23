## code to process Ig GWAS summary statistics so we can run through CHICGP pipeline

library(data.table)

DATA.DIR<-'/home/ob219/scratch/Ig_gwas/raw'
PMI.OUT.DIR<-'/home/ob219/scratch/Ig_gwas/pmi'
GENE.SCORE.OUT.DIR<-'/home/ob219/scratch/Ig_gwas/gene_score'
SH.FILE<-'/home/ob219/scratch/Ig_gwas/computeGeneScores.sh'
info.impute.thresh<-0.4
maf.thresh<-0.01
## remove those with HWE p.val < chi2=25
hwe.thresh<-2*pnorm(sqrt(25),lower.tail=FALSE)

fs<-list.files(path=DATA.DIR,pattern="*.txt",full.names=TRUE)
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

Var.data <- function(f, N) {
  1 / (2 * N * f * (1 - f))
}
## compute approx bayes factors and resultant posterior probabilities
## based on the assumption of one causal variant in a region
approx.bf <- function(z,f, N, pi_i=1e-4) {
  sd.prior <- 0.15
  V <- Var.data(f, N)
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

processFile<-function(f){
  DT<-fread(f)
  DT<-subset(DT,hwe>hwe.thresh & maf>maf.thresh & info_impute > info.impute.thresh)
  DT<-addLDBlock(DT,ld.gr)
  DT[,ppi:=approx.bf(beta/se,maf,N),by=ld]
  ## this can then be input into the main pipeline if we output to pmi format
  out.DT<-with(DT,data.table(chr=chr,start=pos-1,end=pos,rsid=rsid,maf=maf,pval=p_value,ppi=ppi,imp.snp.pos=0,imp.r2=1))
  options(scipen=999)
  pmi.o.file<-file.path(PMI.OUT.DIR,paste(sub("\\.txt","",basename(f)),'pmi',sep='.'))
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

sh_cmd<-lapply(fs,processFile)

write(unlist(sh_cmd),file=SH.FILE)
