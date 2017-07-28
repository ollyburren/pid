library(data.table)

#DT<-fread("/home/ob219/scratch/Ig_gwas/LN_IgG_total.txt")
# save in RDS format for future loading and manipulation
#saveRDS(DT,file="/home/ob219/scratch/Ig_gwas/LN_IgG_total.RDS")

DT<-readRDS("/home/ob219/scratch/Ig_gwas/LN_IgG_total.RDS")
# what is n ?

max(DT$N)

## assume looking at variants above 1% MAF.

DT.f<-subset(DT,maf>=0.01 & hwe>1e-4)

## compute genomic inflation

DT.f[,Z:=beta/se]
DT.f[,chi.sq:=Z^2]

lambda<-mean(DT.f$chi.sq)

## to adjust

DT.f[,p.adj:=2 * pnorm(sqrt(chi.sq/lambda),lower.tail=FALSE)]

# next plot qqplot

QQ=function(x,l=0.99,add=FALSE,minx=TRUE,...) {
  if (max(abs(x),na.rm=T)<1.1) x=-log10(x)
  n=length(x); q=-log10((n:1)/(n+1)); x=sort(x)
  n1=round(l*n)
  if (minx) {
    if (add) points(c(0,q[n1:n]),c(0,x[n1:n]),...) else plot(c(0,q[n1:n]),c(0,x[n1:n]),...)
    lines(c(0,q[n1]),c(0,x[n1]),...)
  } else if (add) points(q,x,...) else plot(q,x,...)
}

QQ(-log10(DT.f$p.adj))

## what happens if we remove MHC

## code for plotting a Manhattan plot in ggplot

library(ggplot2)

DT.f<-DT.f[order(DT.f$chr,DT.f$pos),]
tmp<-split(DT.f$pos,DT.f$chr)


cs<-c(0,head(cumsum(as.numeric(sapply(tmp,max))),-1))
for(i in seq_along(tmp)){
  tmp[[i]]<-tmp[[i]] + cs[i]
}

DT.f$gpos<-do.call('c',tmp)

## for manhatten plotting we don't care about p-values under a certain threshold
THRESHOLD<-1e-2
ggplot(DT.f[DT.f$p.adj<THRESHOLD],aes(x=gpos,y=-log10(p.adj),color=chr%%2==0)) + geom_point(size=0.1) + theme_bw() +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
        geom_hline(yintercept = c(-log10(5e-8),-log10(1e-5)),color=c('red','blue')) + coord_cartesian(ylim=c(0,25))

## investigate the hits

hits<-subset(DT.f,p.adj<5e-8)

library(GenomicRanges)

hits.gr<-with(hits,GRanges(seqnames=Rle(chr),ranges=IRanges(start=pos,end=pos),p=p.adj,name=rsid))
## load in recomb blocks for rough idea
ld<-fread("/home/ob219/scratch/DATA/JAVIERRE_GWAS/support/0.1cM_regions.b37.bed")
ld[,V2:=V2+1]
ld.gr<-with(ld,GRanges(seqnames=Rle(V1),ranges=IRanges(start=V2,end=V3),id=V4))
mergeByOverlaps(ld.gr,hits.gr)
hits<-mergeByOverlaps(ld.gr,hits.gr)
top.hits<-do.call('rbind',lapply(split(hits,hits$id),function(h){
  h[which.min(h$p),]
}))

# check differential missingness



## what we want is to show in IGV and then we can send around
for.igv<-DT.f[DT.f$rsid %in% hits$name,.(chr,pos,rsid,p.adj)]
setnames(for.igv,c('CHR','POS','SNP','P'))
write.table(for.IGV,file='/home/ob219/scratch/Ig_gwas/hits/IgG_total.gwas',sep="\t",row.names=FALSE,quote=FALSE)

## prep and create a track for Jonson et al results - downloaded supp table from https://www.nature.com/ng/journal/v49/n8/extref/ng.3897-S6.xlsx
## save IgG results as csv file

dat<-fread("IgG_jonsson.csv")
dat<-dat[,c(1:2,14),with=FALSE]
setnames(dat,c('SNP','POS','P'))
dat<-cbind(dat,do.call('cbind',(tstrsplit(dat$POS, ":", fixed=TRUE))))
dat$POS<-NULL
setnames(dat,c('SNP','P','CHR','POS'))
setcolorder(dat,c('CHR','POS','SNP','P'))
dat$POS<-as.numeric(dat$POS)

## move to hg19
library(rtracklayer)
chain.file.path<-'/Users/oliver/Downloads/hg38ToHg19.over.chain'
c<-import.chain(chain.file.path)
dat.gr<-with(dat,GRanges(seqnames=Rle(CHR),ranges=IRanges(start=POS,end=POS),SNP=SNP,P=P))
example.37.gr<-unlist(liftOver(dat.gr,c))
names(mcols(example.36.gr))<-'snp.name'
dat.37<-data.table(as.data.frame(example.37.gr))
dat.37<-dat.37[,.(seqnames,start,SNP,P)]
setnames(dat.37,c('CHR','POS','SNP','P'))
dat.37$CHR<-sub("^chr","",dat.37$CHR)

## shell script to process Rafal's slides.

library(data.table)
library(magrittr)
analysis<-'SNPTEST'
an.reg<-sprintf("/%s/",analysis)
fs<-list.files(path="/scratch/ru222/old/data_20170104-B/",recursive=TRUE,full.names=TRUE)
fs<-fs[grep(an.reg,fs,fixed=TRUE)]
## remove the txt files
fs<-fs[grep('txt$',fs,invert=TRUE)]
names(fs)<-sub(".*(chr[0-9]+).*","\\1",fs)
pid.res<-lapply(fs,function(f){
  message(sprintf("Processing %s",f))
  chr<-sub(".*chr([0-9]+).*","\\1",f)
  cmd<-sprintf("grep -v '^#' %s | cut -d' ' -f4,31,32,33,35,45,47,48",f)
  DT<-fread(cmd)
  DT[,chr:=chr]
  DT[,chi.sq:=(frequentist_add_beta_1/frequentist_add_se_1)^2]
})

pid.res<-rbindlist(pid.res)

## filter so in line with other gwas
pid.res.f<-subset(pid.res,controls_maf>=0.01 & cohort_1_hwe>1e-4 & missing_data_proportion <0.02)

## remove HWE

## adjust for population inflation

lambda.pid<-mean(pid.res.f$chi.sq)
## adjust p values accordingly
pid.res.f[,p.adj:=2 * pnorm(sqrt(chi.sq/lambda),lower.tail=FALSE)]

## combine two data sets
DT.f[,id:=paste(chr,pos,sep=':')]
setkey(DT.f,id)
pid.res.f[,chr:=sub("chr","",chr,fixed=TRUE)]
pid.res.f[,id:=paste(chr,position,sep=':')]
setkey(pid.res.f,id)

out<-DT.f[pid.res.f]

## remove MHC 6:31343101-32804570 and chr14 14:105955669-106397775 region

rm.idx<-with(out,which((between(pos,31e6,35e6) & chr==6) | between(pos,105955669,106397775) & chr==14))

## quick look

q<-out[-rm.idx,]

CQ=function (xi, xj, cuts = 10^-(1:5), l = 0.99, minx = TRUE, ...)
{
    if (!(length(xi) == length(xj))) 
        stop("Parameters xi and xj must be vectors of the same length")
    w = which(is.finite(xi + xj))
    xi = xi[w]
    xj = xj[w]
    if (max(abs(xi), na.rm = T) < 1.1)
        xi = -log10(xi)
    if (max(abs(xj), na.rm = T) < 1.1)
        xj = -log10(xj)
    QQ(xi, l = l, minx = minx, col = "black", ...)
    for (i in 1:length(cuts)) {
        sub = which(xj > -log10(cuts[i]))
        if (length(sub) * (1 - l) > 100)
            QQ(xi[sub], l = l, minx = minx, col = i + 1, add = T,
                ...)
        if ((length(sub) > 10) & (length(sub) * (1 - l) < 100))
            QQ(xi[sub], minx = F, col = i + 1, add = T, ...)
    }
    abline(0, 1, col = "red", lty = 2)
}

CQ(q$i.p.adj,q$p.adj)
