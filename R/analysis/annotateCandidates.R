library(data.table)
library(GenomicRanges)

genes<-unique(c('ENSG00000134352',
'ENSG00000056558',
'ENSG00000134460',
'ENSG00000073050',
'ENSG00000101096',
'ENSG00000141564',
'ENSG00000116824',
'ENSG00000116824',
'ENSG00000116824',
'ENSG00000117322',
'ENSG00000198589',
'ENSG00000198589',
'ENSG00000198589'))

individuals<-unique(c('F010564',
'F010256',
'F011392',
'F009298',
'F009963',
'F009873',
'F009297',
'F009297',
'F009297',
'F013197',
'F000929',
'F000929',
'F000929'))


DATA.DIR<-'/home/ob219/scratch/pid/ANNOTATIONS/'

library(magrittr)

prom<-readRDS(file.path(DATA.DIR,'promoters.RDS')) %>% subset(.,ensg %in% genes)

## exon regions
exon<-readRDS(file.path(DATA.DIR,'exons.RDS')) %>% subset(.,ensg %in% genes)

## pchic superenhancer (for time being)
pse<-readRDS(file.path(DATA.DIR,'pchic_hnisz','se.gr.RDS')) %>% subset(.,ensg %in% genes)

genes<-readRDS(file.path(DATA.DIR,'genes.RDS')) %>% subset(.,ensg %in% genes)

manta.files<-list.files(path='/home/ob219/scratch/pid/MANTA/',pattern="*.RDS",full.names=TRUE)
mevents<-lapply(manta.files,readRDS)
names(mevents)<-gsub("pid\\_([^_]+).*","\\1",basename(manta.files))
mevents<-lapply(mevents,function(m){
  subset(m,bridge_id %in% individuals)
})
## drop bnd for time being
mevents<-mevents[names(mevents) != 'bnd']

computeFeatOl<-function(e.gr,f.gr){
  out<-logical(length=length(e.gr))
  ol<-unique(as.matrix(findOverlaps(e.gr,f.gr))[,1])
  out[ol]<-TRUE
  out
}

fl<-list(prom=prom,exon=exon,pse=pse)

annoMatrix<-function(e){
  tmp<-do.call('cbind',lapply(fl,function(f.gr){
    computeFeatOl(e,f.gr)
  }))
  ## get a list of all those events that overlap at least one features
  keep<-which(rowSums(tmp)!=0)
  info<-cbind(mcols(e),DataFrame(tmp))[keep,]
  e<-e[keep,]
  mcols(e)<-info
  e
}

res<-lapply(mevents,annoMatrix)

## for these want feature genes (i.e. a denormalised list of CNV's and genes and reasons)

mergeFeatures<-function(f,flist=fl){
  message("Merging features")
  mf<-mcols(f)
  rbindlist(lapply(names(fl),function(tname){
    message(sprintf("Merging %s",tname))
    p<-f[mf[[tname]],]
    m<-data.table(as.data.frame(mergeByOverlaps(p,fl[[tname]])))
    m<-m[,.(p.seqnames,p.start,p.end,p.gt,p.uid,p.bridge_id,ensg,name)]
    m$ftype<-tname
    m
  }))
}

all<-lapply(res,function(x) unique(mergeFeatures(x)))

snps<-rbindlist(lapply(list.files(path='/home/ob219/scratch/pid/VCF/parse_genes/RDS',full.names=TRUE),readRDS))

## merge the events with related SNPs to get a really large list
pid.genes<-fread("/home/ob219/scratch/pid/all.pid.genes.csv")
poss.diag<-c('F000891','F000891','F000903','F000913','F000923','F000924','F000932','F000940','F000943','F000956','F000958','F001339','F001341','F001342','F001356','F003358','F003360','F003374','F003380','F003449','F003452','F003458','F003462','F003464','F003464','F004692','F005170','F005170','F005216','F005217','F005220','F005474','F005474','F005482','F005484','F005533','F005987','F005989','F005990','F006326','F006352','F006352','F006364','F006365','F006367','F006422','F006422','F006467','F006478','F006484','F006486','F006492','F007088','F007163','F007163','F007164','F007164','F007179','F008103','F008815','F008870','F008871','F008876','F008877','F008877','F008988','F009202','F009220','F009308','F009726','F009988','F010077','F010232','F010257','F010336','F010350','F010373','F010381','F010382','F010434','F010560','F010562','F010599','F010732','F010738','F010773','F010773','F011022','F011419','F012317','F012524','F013187','F013188','F013235','F014738','F014738','W000010','W000011','W000232','W000272','W000288','F000908','F000938','F001274','F005527','F005530','F005532','F006340','F006496','F008512','F008970','F009790','F009795','F011386','F013197','F005212','F005212','F009729','F010048','F010052','F010052','F010080','F010731','F011393')
def.diag<-c('F000891','F000903','F000913','F000923','F000924','F000932','F000940','F000943','F000956','F000958','F001339','F001341','F001342','F001356','F003358','F003360','F003374','F003380','F003449','F003452','F003458','F003462','F003464','F004692','F005170','F005216','F005217','F005220','F005474','F005482','F005484','F005533','F005987','F005989','F005990','F006326','F006352','F006364','F006365','F006367','F006422','F006467','F006478','F006484','F006486','F006492','F007088','F007163','F007164','F007179','F008103','F008815','F008870','F008871','F008876','F008877','F008988','F009202','F009220','F009308','F009726','F009988','F010077','F010232','F010257','F010336','F010350','F010373','F010381','F010382','F010434','F010560','F010562','F010599','F010732','F010738','F010773','F011022','F011419','F012317','F012524','F013187','F013188','F013235','F014738','W000010','W000011','W000232','W000272','W000288')
ws<-lapply(names(all),function(n){
    x<-all[[n]]
    tmp<-unique(merge(x,snps,by.x=c('ensg','p.bridge_id'),by.y=c('Gene','ind'),allow.cartesian=TRUE))
    tmp<-subset(tmp,BIOTYPE=='protein_coding' & gt %in% c(2,3) & OPR>0.8)
    tmp<-unique(tmp[,.(ensg,name,ftype,p.bridge_id,p.seqnames,p.start,p.end,p.gt,gt,GNOMAD_AF,WGS10K_AF,Existing_variation,Consequence,CADD_PHRED)])
    tmp$pid<-tmp$ensg %in% pid.genes$ensembl_gene_id
    tmp$sv.gt<-0
    tmp[tmp$p.gt %in% c('1/0','0/1')]$sv.gt<-2
    tmp[tmp$p.gt=='1/1']$sv.gt<-3
    tmp[tmp$p.gt=='0/0']$sv.gt<-1
    tmp$p.gt<-NULL
    tmp$sv.type<-n
    setnames(tmp,c('ensg','gene','feat.type','individual','chr','sv.start','sv.end','snp.gt','gnomad_af','wgs10k_af','exist.var','consequence','cadd','pid','sv.gt','sv.type'))
    setcolorder(tmp,c('sv.type','ensg','gene','feat.type','individual','chr','sv.start','sv.end','sv.gt','snp.gt','gnomad_af','wgs10k_af','exist.var','consequence','cadd','pid'))
    tmp$poss.diag<-tmp$individual %in% poss.diag
    tmp$def.diag<-tmp$individual %in% def.diag
    ##fix numeric columns
    for(n in c('snp.gt','gnomad_af','wgs10k_af','cadd')){
      tmp[[n]]<-as.numeric(tmp[[n]])
    }
    tmp
})

## filters

## take a look at frameshift

getByCons<-function(DT,term){
  DT[grep(term,DT$consequence),]
}

filtByAF<-function(DT,GAF,WAF){
  DT[(is.na(DT$gnomad_af) | as.numeric(DT$gnomad_af)<GAF) & DT$wgs10k_af<WAF,]
}

filtByCADD<-function(DT,score){
  DT[as.numeric(DT$cadd)>score,]
}

filtByType<-function(DT,type){
  DT[DT$feat.type==type,]
}

filtByPID<-function(DT){
  DT[DT$pid,]
}


all.rare<-lapply(ws,filtByAF,GAF=1e-3,WAF=1e-3)
all.hi.cadd.rare<-lapply(all.rare,filtByCADD,score=20)

by.pid<-rbindlist(lapply(all.hi.cadd.rare,filtByPID))

final.list<-rbindlist(all.hi.cadd.rare)

library(xlsx)
write.xlsx(final.list, file="/home/ob219/scratch/pid/RESULTS/possible_CHET.xlsx", sheetName="all")

## add more information for pse

pses<-subset(final.list,feat.type=='pse')

pses.gr<-with(pses,GRanges(seqnames=Rle(chr),ranges=IRanges(start=sv.start,end=sv.end),individual=individual))

pse.anno<-as.data.frame(mergeByOverlaps(pses.gr,pse))

write.xlsx(pse.anno, file="/home/ob219/scratch/pid/RESULTS/possible_CHET.xlsx", sheetName="pse_cell_types")

stop()


## process SV homs


homs<-rbindlist(lapply(names(all),function(n){
  tmp<-subset(all[[n]],p.gt=='1/1')
  tmp$sv.gt<-0
  tmp[tmp$p.gt %in% c('1/0','0/1')]$sv.gt<-2
  tmp[tmp$p.gt=='1/1']$sv.gt<-3
  tmp[tmp$p.gt=='0/0']$sv.gt<-1
  tmp$p.gt<-NULL
  tmp$sv.type<-n
  tmp$poss.diag<-tmp$p.bridge_id %in% poss.diag
  tmp$def.diag<-tmp$p.bridge_id %in% def.diag
  tmp$pid<-tmp$ensg %in% pid.genes$ensembl_gene_id
  tmp<-tmp[,.(sv.type,ensg,name,ftype,p.bridge_id,p.seqnames,sv.gt,pid,poss.diag,def.diag)]
  setnames(tmp,c('sv.type','ensg','gene','feat.type','individual','chr','sv.gt','pid','poss.diag','def.diag'))
  tmp$chr<-as.numeric(tmp$chr)
  tmp
}))

write.xlsx(homs, file="/home/ob219/scratch/pid/RESULTS/hom_sv.xlsx", sheetName="all")

stop()

hets<-lapply(all,function(x) subset(x,p.gt!='1/1'))
hom<-lapply(all,function(x) subset(x,p.gt=='1/1'))

## there are lots of SV events overlaping PCHIC superenhancers - we can
## attempt to prioritise by looking at those genes that are affected by exon / prom effects in some individuals
## and pse in others as the different evidence types might be nice.
## of course by nature this excludes events occuring in a single individual.
hom.mix.cause<-lapply(hom,function(test){
  #split by gene
  if(nrow(test)==0)
    return(NA)
  t2<-split(test,test$ensg)
  ## there a lots of PSE what might be interesting is prom / exon combined with a pse
  interest<-sapply(t2,function(x){
    sum(x$ftype!='pse') >0 & sum(x$ftype=='pse') > 0
  })
  unique(t2[interest])
})

## inv - all on F000929 - large inversion event on chromosome X is M or F
## del - none - get a list of genes affecting exon and prom
## dup - no homs in the first place
## ins - 9 genes - individuals per gene - need to get genes and number of individuals
interesting.hom.ins<-lapply(hom.mix.cause[['ins']],function(x){
  dups<-x[duplicated(x$p.bridge_id),]$p.bridge_id
  tmp<-x[!x$p.bridge_id %in% dups,]
  if(sum(tmp$ftype!='pse') >0){
    return(tmp)
  }
  return(NA)
})

## TODO - look and see if any of these also associated with LoF variaition in target gene.
## need a pipeline for this.

## look at heterozygotes.
## first look at deletions and see if we can find arpc1b example by includinging genotypes.


dels<-hets[['del']]

## to find ARPC1B we first get all genes that have only a single case when pse is included.
sdel<-dels[,list(scount=length(unique(p.bridge_id)),
  inds=paste(unique(p.bridge_id),sep=',',collapse=','),
  event.type=paste(unique(ftype),sep=',',collapse=',')),by=ensg]
## events that occur in one individual for a gene
ug<-subset(sdel,scount==1)

## ok filter this list by those that have lof in the same gene

snps<-readRDS("/home/ob219/scratch/pid/ANNOTATIONS/lof_snps.RDS")
## make the above by individual and gene.
snps.f<-subset(snps,gene %in% ug$ensg)
snps.f$id<-1:nrow(snps.f)
gt<-strsplit(snps.f$gt,',')

inds<-rbindlist(lapply(seq_along(gt),function(i){
  mgt<-gt[[i]]
  tmp<-data.table(do.call('rbind',strsplit(mgt,':')))
  setnames(tmp,c('bridge_id','geno'))
  tmp$id=i
  tmp
}))

m<-merge(snps.f,inds,by.x='id',by.y='id')
m$gt<-NULL

## merge snps with SV events

fm<-merge(ug,m,by.x=c('ensg','inds'),by.y=c('gene','bridge_id'))

## just get rare snps

fm.rare<-subset(fm,wgs10k<1e-3)

## annotate these with pid status genes

pid.genes<-fread("/home/ob219/scratch/pid/all.pid.genes.csv")
fm.rare$pid<-fm.rare$ensg %in% pid.genes$ensembl_gene_id

## reveals 4 individuals with possible CHETs in known pid genes - need to verify. - they are all het at possible
## LOF mutations and the genotypes look rare - need to visualise to show everyone.

saveRDS(subset(fm.rare, pid),"/home/ob219/scratch/pid/RESULTS/deletion_rare_chet_known_pid.RDS")

## looked at these in validate_pid_chets.R - deletions look good but SNPs are not that convincing
## possible SE deletion in conjunction with missense in DCLRE1B.
## possible prom,exon deletion in conjunction with UTR variant in STAT1 (possible dominiant phenotype otherwise ?)
## need to check if anyone else has a deletion affecting the same gene.

saveRDS(subset(fm.rare, !pid),"/home/ob219/scratch/pid/RESULTS/deletion_rare_chet_known_NOT_pid.RDS")
##
