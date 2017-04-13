library(biomaRt)
library(data.table)

data.dir<-'/scratch/ob219/pid/'
## add in ensembl id's
pid.genes<-fread(file.path(data.dir,'Tier1-IUIS_gene_transcripts_sort.txt'))
setnames(pid.genes,c('chr','start','end','gene','anno','ense','size'))
pid.genes.filt<-subset(pid.genes,ense != '.')

mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
ensembl.gene <- useDataset("hsapiens_gene_ensembl",mart=mart)

gene.details<-getBM(
    filters= c("ensembl_exon_id"),
    attributes= c('ensembl_gene_id','ensembl_exon_id','external_gene_name'),
    values= list(ensembl_gene_id=unique(pid.genes.filt$ense)),
    mart= ensembl.gene)

## load in the full list to check that the things that we expect are found
pid.genes.name<-scan(file.path(data.dir,'PID-Tier1-IUIS-genes.txt'),character())
idx<-which(pid.genes.name %in% gene.details$external_gene_name)
miss.pid.genes<-pid.genes.name[-idx]
## are these actually in transcript file ?
sum(miss.pid.genes %in% pid.genes$gene)
# no can we get from biomart by external_gene_name ?
gene.details2<-getBM(
    filters= c("external_gene_name"),
    attributes= c('ensembl_gene_id','ensembl_exon_id','external_gene_name'),
    values= list(ensembl_gene_id=miss.pid.genes),
    mart= ensembl.gene)

all<-rbind(gene.details,gene.details2)
idx<-which(pid.genes.name %in% all$external_gene_name)
## annotate these by hand
all<-unique(data.table(all[,c('ensembl_gene_id','external_gene_name')],key='ensembl_gene_id'))
all$pid.name<-all$external_gene_name
manual<-data.table(ensembl_gene_id=c('ENSG00000093072','ENSG00000177084',NA),external_gene_name=c('CECR1','POLE',NA),pid.name=c('ADA2','POLE1','MHS6'))
all.pid<-rbind(all,manual)
write.table(all.pid,file=file.path(data.dir,'all.pid.genes.csv'),sep=',',quote=FALSE,row.names=FALSE)
## are all PID genes found ab negative
ab.def<-fread(file.path(data.dir,'ab_def_genes.csv'),select=1:4)
ab.def$ENSG %in% all.pid$ensembl_gene_id
## next read in PCHiC Dataset so that we can filter and get lists of PIRs that interact with these genes. Should we also think
## about promoter regions ?
chic<-fread(file.path(data.dir,"merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab"))
chic.pid<-subset(chic, ensg %in% all.pid$ensembl_gene_id)
## which genes have no PIR in our cell types
pid.genes.interaction<-all.pid[all.pid$ensembl_gene_id %in% chic.pid$ensg,]$ensembl_gene_id
## whats the missing genes for which we have no results - perhaps a missing bait ?
idx<-which(all.pid$ensembl_gene_id %in% pid.genes.interaction)
pid.genes.no.interaction<-all.pid[-idx,]
save(chic.pid,file=file.path(data.dir,'pid.chic.RData'))

## what about all genes rather than just PID ? Can just use chic<-fread(file.path(data.dir,"merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab"))
