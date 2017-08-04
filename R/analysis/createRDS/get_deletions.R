library(data.table)
library(GenomicRanges)

DATA.DIR<-'/scratch/WGS10K/data/release/latest/merged-sv/all/'
OUT.DIR<-'/home/ob219/scratch/pid/MANTA/'
FREQ<-0.03

## get list of ID's to include in analysis1

pid.samples<-scan('/home/ob219/scratch/pid/support/PID_848_index_cases_no_vasculitis_3_8_2017.txt',character())

## get list of those with diagnosis

diag<-fread("/home/ob219/scratch/pid/support/BRIDGE-PID_DIAGNOSIS_3_8_2017.csv")
## remove those with diag of full - think that this has already been done.
pid.samples<-setdiff(pid.samples,subset(diag,variant.contribution.phenotype=='full')$WGS.ID)

## get MANTA calls


getMANTA<-function(f,samples=pid.samples){
  manta<-fread(f)
  manta[,uid:=1:nrow(manta)]
  manta<-manta[manta$BRIDGE_ID %in% samples,.(BRIDGE_ID,MANTA_CHROM,MANTA_START,MANTA_END,MANTA_GT,MANTA_ALT,uid)]
  with(manta,GRanges(seqnames=Rle(MANTA_CHROM),ranges=IRanges(start=MANTA_START,end=MANTA_END),gt=MANTA_GT,uid=uid,bridge_id=BRIDGE_ID,alt=MANTA_ALT))
}

manta.files<-list.files(path=DATA.DIR,pattern="mantacalls*",full.names=TRUE)
## get just the correct frequency and remove exons etc
manta.files<-manta.files[intersect(grep(FREQ,manta.files),grep("exons",manta.files,invert=TRUE))]

for(f in manta.files){
  message(sprintf("Processing %s",f))
  of<-paste('pid',sub("mantacalls_allsamples\\_([^\\.]+)\\..*","\\1",basename(f)),FREQ,sep='_')
  m<-getMANTA(f)
  saveRDS(m,file=file.path(OUT.DIR,paste(of,'RDS',sep='.')))
}
