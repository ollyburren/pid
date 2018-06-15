## test code to obtain pvalues for enrichment if found.

library(data.table)
library(rtracklayer)
manifest.file <- '/home/ob219/scratch/pid/support/CondGWAS.tsv'
m.DT <- fread(manifest.file)
gwas.dir <- '/home/ob219/rds/hpc-work/DATA/JAVIERRE_GWAS/gwas/'

## process Ab deficient
if(!file.exists('/home/ob219/rds/hpc-work/DATA/JAVIERRE_GWAS/gwas/AB_DEF_AI_COV_2018.bed.gz')){
  ab.def.dir <- '/rds/user/hl454/hpc-work/GWAS/plink_case_control/AbDeficiency2/'
  abdf <- fread('/scratch/WGS10K/analysis/pid/AbDeficiency_covAI.chr1-22.assoc.logistic.plink')

#  abdf <- list.files(path=ab.def.dir,pattern='*.logistic',full.names=TRUE)

#  abdf <- rbindlist(lapply(abdf,function(f){
#    message(f)
#    DT <- fread(f)[,.(CHR,BP-1,BP,SNP,P)]
#    setnames(DT,c('chr','start','end','id','p.val'))
#  }))
  abdf<-abdf[,.(CHR,BP-1,BP,SNP,P)]
  setnames(abdf,c('chr','start','end','id','p.val'))
  abdf <- abdf[order(chr,start),]
  options(scipen=999)
  write.table(abdf,file="/home/ob219/rds/hpc-work/DATA/JAVIERRE_GWAS/gwas/AB_DEF_AI_COV_2018.bed",col.names=FALSE,row.names=FALSE,quote=FALSE)
  options(scipen=0)
}

## process to get just variants at 1% from a  reference data set such as 1K genomes.
#maf.file <- '/home/ob219/rds/hpc-work/DATA/1kgenome/SNPMATRIX/summary_1_per/all.chr.RDS'
#if(!file.exists(maf.file)){
#  library(data.table)
#  files <- list.files(path='/home/ob219/rds/hpc-work/DATA/1kgenome/SNPMATRIX/summary/',pattern='.*',full.names=TRUE)
#  all.chr <- lapply(files,function(f){
#    message(f)
#    LA <- get(load(f))[MAF>0.01,][,chr:=gsub("([^\\.]+)\\..*","\\1",basename(f))]
#    LA <- LA[,.(uid,MAF)]
#    LA[,MAF:=round(MAF,digits=3)]
#    LA[,uid:=as.character(uid)]
#    LA
#    })
#  all.1k <- rbindlist(all.chr)
#  setkey(all.1k,uid)
#  saveRDS(all.1k,file=maf.file)
#}else{
#  all.1k <- readRDS(maf.file)
#}
MAF.thresh <- 0.01 #1%
(load("/home/ob219/scratch/as_basis/1KG_support/all_EUR_support.RData"))
options(scipen=999)
all.eur[,pid:=paste(chr,position,sep=':')]
options(scipen=0)
all.eur[,maf:=a2.f]
all.eur[maf>0.5,maf:=1-maf]
all.eur <- all.eur[,maf:=round(maf,digits=3)]
all.eur <- all.eur[maf>MAF.thresh,]
setkey(all.eur,pid)


ld.gr<-import.bed("/home/ob219/scratch/DATA/JAVIERRE_GWAS/support/0.1cM_regions.b37.bed")

addLDBlock<-function(DT,ld.gr){
  dt.gr<-with(DT,GRanges(seqnames=Rle(chr),ranges=IRanges(start=position,width=1L),idx=1:nrow(DT)))
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


out.dir <- '/home/ob219/scratch/pid/GWAS_AI_COVARIATE'

for(i in 1:nrow(m.DT)){
  fname <- sprintf("%s.RDS",m.DT$label[i])
  if(file.exists(file.path(out.dir,fname)))
    next
  message(sprintf("Processing %s",m.DT$label[i]))
  cmd <- sprintf("zcat %s",file.path(gwas.dir,m.DT$filename[i]))
  DT.f <- fread(cmd)
  setnames(DT.f,c('chr','start','end','id','p.val'))
  DT.f <- DT.f[!is.na(p.val),]
  options(scipen=999)
  DT.f[,id:=paste(chr,end,sep=':')]
  options(scipen=0)
  setkey(DT.f,id)
  DT.maf <- all.eur[DT.f][!is.na(maf),.(name,chr,position,pid,maf,p.val)]
  rm.idx<-with(DT.maf,which((between(position,25e6,45e6) & chr==6)))
  DT.maf<-DT.maf[-rm.idx,]
  #ret[,c('trait','cases','controls':=list(m.DT$label[i],m.DT$cases[i],m.DT$comtrols[i])]
  DT.maf<-addLDBlock(DT.maf,ld.gr)
  ## need to add MAF
  if(is.na(m.DT$controls[i])){
    DT.maf[,ppi:=approx.bf.p(p.val,maf,m.DT$cases[i],0,1e-4,'QUANT'),by=ld]
  }else{
    N<- m.DT$cases[i] + m.DT$controls[i]
    DT.maf[,ppi:=approx.bf.p(p.val,maf,N,m.DT$cases[i]/N,1e-4),by=ld]
  }

  saveRDS(DT.maf,file=file.path(out.dir,fname))
  message(sprintf("Saved to %s",file.path(out.dir,fname)))
  summ <- DT.maf[,list(total=sum(ppi)),by=ld]
  fname <- sprintf("ld_summary_%s.RDS",m.DT$label[i])
  saveRDS(summ,file=file.path(out.dir,fname))
  message(sprintf("Saved to %s",file.path(out.dir,fname)))
}


## work out enrichment or not

bfiles <- list.files(path=out.dir,pattern="ld_summary.*",full.names=TRUE)
all.traits <- lapply(bfiles,readRDS)
names(all.traits)<-gsub("ld_summary_([^.]+)\\..*","\\1",basename(bfiles))
library(wgsea)
pmi.thresh <- 0.8
n.perms <- 1e3
AbDef<- all.traits$ABDEF
all.W <- lapply(seq_along(all.traits),function(i){
  trait <- names(all.traits)[i]
  message(trait)
  DT <- all.traits[[i]]
  ## regions
  R <- DT[total>pmi.thresh,]$ld
  idx <- which(AbDef$ld %in% R)
  W <- wilcox.test(AbDef$total[idx],AbDef$total[-idx])$statistic
  Wstar <- numeric(length=n.perms)
  for(j in 1:n.perms){
    message(j)
    idx.sample <- sample(1:nrow(AbDef),length(R),replace=FALSE)
    Wstar[j] <- wilcox.test(AbDef$total[idx.sample],AbDef$total[-idx.sample])$statistic
  }
  enrich<-Z.value(W=W,Wstar=Wstar,n.in=length(idx),n.out=length(AbDef$total[-idx]))
  data.table(trait=trait,Z=enrich$Z.empirical$statistic,P=enrich$Z.empirical$p.value)
})

all.W <- rbindlist(all.W)
all.W <- all.W[order(all.W$Z),]
all.W <- all.W[!trait %in% c('ABDEF_AI_COV','ABDEF')]
all.W[trait %in% c('RA','SLE','T1D','UC','CD'),fill:='Autoimmune/Autoinflammatory']
all.W[!trait %in% c('RA','SLE','T1D','UC','CD'),fill:='Other']
all.W$trait<-factor(all.W$trait,levels=all.W$trait)

saveRDS(all.W,file="~/tmp/abdef.RDS")

sig.thresh <- qnorm(0.05/length(unique(all.W$trait)),lower.tail=FALSE)

library(ggplot2)
pdf(file="~/tmp/abdef_enrich.pdf")
ggplot(all.W,aes(x=trait,y=Z,fill=fill)) + geom_bar(stat="identity") + theme_bw() +
scale_fill_manual(name='Category',labels=c("Autoimmune/\nAutoinflammatory","Other"),values=c('firebrick','dodgerblue')) +
xlab("Trait") + ylab("Empirical Wilcoxon Z score") + geom_hline(yintercept=sig.thresh,color='red',lty=2)
dev.off()


bfiles <- list.files(path=out.dir,pattern="ld_summary.*",full.names=TRUE)
all.traits <- lapply(bfiles,readRDS)
names(all.traits)<-gsub("ld_summary_([^.]+)\\..*","\\1",basename(bfiles))
library(wgsea)
pmi.thresh <- 0.8
n.perms <- 1e3
AbDef<- all.traits$ABDEF_AI_COV
all.W <- lapply(seq_along(all.traits),function(i){
  trait <- names(all.traits)[i]
  message(trait)
  DT <- all.traits[[i]]
  ## regions
  R <- DT[total>pmi.thresh,]$ld
  idx <- which(AbDef$ld %in% R)
  W <- wilcox.test(AbDef$total[idx],AbDef$total[-idx])$statistic
  Wstar <- numeric(length=n.perms)
  for(j in 1:n.perms){
    message(j)
    idx.sample <- sample(1:nrow(AbDef),length(R),replace=FALSE)
    Wstar[j] <- wilcox.test(AbDef$total[idx.sample],AbDef$total[-idx.sample])$statistic
  }
  enrich<-Z.value(W=W,Wstar=Wstar,n.in=length(idx),n.out=length(AbDef$total[-idx]))
  data.table(trait=trait,Z=enrich$Z.empirical$statistic,P=enrich$Z.empirical$p.value)
})

all.W <- rbindlist(all.W)
all.W <- all.W[order(all.W$Z),]
all.W <- all.W[!trait %in% c('ABDEF_AI_COV','ABDEF')]
all.W[trait %in% c('RA','SLE','T1D','UC','CD'),fill:='Autoimmune/Autoinflammatory']
all.W[!trait %in% c('RA','SLE','T1D','UC','CD'),fill:='Other']
all.W$trait<-factor(all.W$trait,levels=all.W$trait)

saveRDS(all.W,file="~/tmp/abdef_ai_cov.RDS")

## what is bf threshold on Z score

sig.thresh <- qnorm(0.05/length(unique(all.W$trait)),lower.tail=FALSE)

library(ggplot2)
pdf(file="~/tmp/abdef_enrich_cov_ai.pdf")
ggplot(all.W,aes(x=trait,y=Z,fill=fill)) + geom_bar(stat="identity") + theme_bw() +
scale_fill_manual(name='Category',labels=c("Autoimmune/\nAutoinflammatory","Other"),values=c('firebrick','dodgerblue')) +
xlab("Trait") + ylab("Empirical Wilcoxon Z score") + geom_hline(yintercept=sig.thresh,color='red',lty=2)
dev.off()



## TODO TAKE KNOWN PID GENES AND THEN WORK OUT IF THERE IS ENRICHMENT OF GWAS SIGNALS IN REGIONS WITH KNOWN PID GENS IN.
