library(data.table)
## take a look at summary table
stable<-fread("/scratch/WGS10K/data/release/latest/merged-sv/summary_table/summary_table_svs_20161012.txt")
## 8 PID samples are 'flagged' perhaps exclude from analysis

## what about strict deletion calls ?
sdir<-'/scratch/WGS10K/data/WGS10K/us/SV-processing-V2all-upgrades/output/output_files/20170104/dels_strict/'
sfiles.bed<-list.files(sdir,pattern='*.bed',full.names=TRUE)
sfiles.txt<-list.files(sdir,pattern='*.txt',full.names=TRUE)

strict.txt<-lapply(sfiles.txt,fread)
names(strict.txt)<-basename(sfiles.txt)
## same format ?
sapply(strict.txt,function(x) length(names(x)))
## matrix files are different assume these are samples on one axis and exons on other
## also deletion.strict.allsamples.dedup.20161012.txt - missing annotations back to genes 
## how many in each file 
sapply(strict.txt,function(x) nrow(x))
## what do the thresholds refer to  - from Olga seems as if these are frequency thresholds about how often these are seen across the cohort. 
all<-strict.txt[['deletions.strict.allsamples.dedup.20161012.txt']]
## what do CANVAS and MANTA FILTERS refer to ? Only make it through if manta and canvas are '.'
## perhaps a good place to start is to look for loss overlap with know PID genes and see how this compares with previous analysis ?
