library(data.table)

source("/home/ob219/git/pid/R/analysis/central_functions.R")

## verify novel CHETs in known genes using raw source files.

candidates<-readRDS("/home/ob219/scratch/pid/RESULTS/deletion_rare_chet_known_pid.RDS")

dels<-fread('/scratch/WGS10K/data/release/latest/merged-sv/all/mantacalls_allsamples_del.ann.filt.0.03.txt')

dels<-subset(dels,BRIDGE_ID %in% candidates$inds & !is.na(ensembl_gene_IDs))

val<-rbindlist(lapply(1:nrow(candidates),function(i){
  DT<-candidates[i,]
  idx<-grep(DT$ensg,dels$ensembl_gene_IDs)
  DT$del.match<-idx
  DT
}))

## cannot verify SE but that makes sense. The others check out.

raw.dels<-dels[val$del.match,]

lapply(1:nrow(val),function(i){
  v<-val[i,]
  print(v)
  r<-sprintf("%s:%s-%s",v$chr,v$position,v$position)
  getPossDamSNPs(r,v$inds,egrep="^#|HGMD_CLASS=DM[^?]|splice|frame|start|stop|missense|utr")
})
