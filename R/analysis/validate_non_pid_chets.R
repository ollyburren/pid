library(data.table)

source("/home/ob219/git/pid/R/analysis/central_functions.R")

## verify novel CHETs in known genes using raw source files.

candidates<-readRDS("/home/ob219/scratch/pid/RESULTS/deletion_rare_chet_not_pid.RDS")

dels<-fread('/scratch/WGS10K/data/release/latest/merged-sv/all/mantacalls_allsamples_del.ann.filt.0.03.txt')

dels<-subset(dels,BRIDGE_ID %in% candidates$inds & !is.na(ensembl_gene_IDs))

val<-rbindlist(lapply(1:nrow(candidates),function(i){
  DT<-candidates[i,]
  idx<-grep(DT$ensg,dels$ensembl_gene_IDs)
  DT$del.match<-idx
  DT
}))

## do some filtering on GNOMAD to remove wacky stuff
val$gnomad_af<-as.numeric(val$gnomad_af)
val[is.na(val$gnomad_af),]$gnomad_af=0
val[,chet_af:=as.numeric(gnomad_af) * 0.03]


## cannot verify SE but that makes sense. The others check out.
raw.dels<-dels[val$del.match,]

funcs<-lapply(1:nrow(val),function(i){
  v<-val[i,]
  print(v)
  r<-sprintf("%s:%s-%s",v$chr,v$position,v$position)
  sn<-getPossDamSNPs(r,v$inds,egrep="^#|HGMD_CLASS=DM[^?]|splice|frame|start|stop|missense")
  cbind(v,sn)
})

## some of these don't pan out due to initial filtering
rm.idx<-which(sapply(funcs,function(j) length(names(j)))!=25)
funcs.f<-rbindlist(funcs[-rm.idx])

val.f<-val[-rm.idx,]
## to prioritise look for single gene events
sing<-val.f[!val.f$ensg %in% val.f[duplicated(val.f$ensg),]$ensg,]



sing.funcs<-rbindlist(lapply(1:nrow(sing),function(i){
  v<-sing[i,]
  print(v)
  r<-sprintf("%s:%s-%s",v$chr,v$position,v$position)
  sn<-getPossDamSNPs(r,v$inds,egrep="^#|HGMD_CLASS=DM[^?]|splice|frame|start|stop|missense")
  cbind(v,sn)
}))

sing.funcs[grep("stop|frameshift",sing.funcs$Consequence),]

top<-subset(sing.funcs,BIOTYPE=='protein_coding' & chet_af<1e-4 & as.numeric(CADD_PHRED)>15)

## yes we found ARPC1B events !!!

subset(top,event.type=='prom,exon')
