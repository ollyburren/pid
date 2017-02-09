library(data.table)
library(GenomicRanges)

chic<-fread(file.path(data.dir,"merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab"))
data.dir<-'/scratch/ob219/pid'
mfiles<-list.files(path=file.path(data.dir,'hnisz'),pattern='*.csv',full.names = TRUE)
## select CD files first off
cd.files<-mfiles[grep('CD',basename(mfiles))]
desc<-vector(mode="character",length=length(cd.files))
desc[grep('CD34',basename(cd.files))]<-'Endothelial_precursors'
## CD4
desc[grep('CD4_Memory',basename(cd.files))]<-'Total_CD4_MF,Total_CD4_Activated,Total_CD4_NonActivated'
desc[grep('BI_CD4p_CD225int_CD127p_Tmem',basename(cd.files))]<-'Total_CD4_MF,Total_CD4_Activated,Total_CD4_NonActivated'
desc[grep('BI_CD4p_CD25-_CD45ROp_Memory',basename(cd.files))]<-'Total_CD4_MF,Total_CD4_Activated,Total_CD4_NonActivated'
desc[grep('CD4_Naive',basename(cd.files))]<-'Total_CD4_MF,Total_CD4_Activated,Total_CD4_NonActivated,Naive_CD4'
desc[grep('BI_CD4p_CD25-_CD45RAp_Naive',basename(cd.files))]<-'Total_CD4_MF,Total_CD4_Activated,Total_CD4_NonActivated,Naive_CD4'
## stimulated T-Cells
desc[grep('BI_CD4p_CD25-_Il17-_PMAstim_Th',basename(cd.files))]<-'Total_CD4_Activated'
desc[grep('BI_CD4p_CD25-_Il17p_PMAstim_Th17',basename(cd.files))]<-'Total_CD4_Activated'
## CD8
desc[grep('BI_CD8_Memory',basename(cd.files))]<-'Total_CD8,Naive_CD8'
desc[grep('CD8_primary',basename(cd.files))]<-'Total_CD8,Naive_CD8'
desc[grep('BI_CD8_Naive',basename(cd.files))]<-'Naive_CD8'
desc[grep('CD14',basename(cd.files))]<-'Macrophages_M0,Macrophages_M1,Macrophages_M2,Monocytes'
desc[grep('CD19',basename(cd.files))]<-'Naive_B,Total_B'
## pro thymocytes
desc[grep('CD3\\.',basename(cd.files))]<-'Naive_B,Total_B,Total_CD4_MF,Total_CD4_Activated,Total_CD4_NonActivated,Naive_CD4,Total_CD8,Naive_CD8'
hnisz<-data.table(f=cd.files,pchic.tissue=desc)
##remove those that don't have a mapping
hnisz<-subset(hnisz,pchic.tissue!='')

pchic.cutoff<-5
pchic<-get(load(file.path(data.dir,'pid.chic.RData')))
mergeSEandPCHiC<-function(f,ti,pchic,se=1){
	DT<-fread(f,skip=1L)
	setnames(DT,c('enh_id','chr','start','end','gene','enh_rank','super','chip.den','read.den'))
	## first off retrieve SE that are found in matched tissues.
	se.gr<-with(DT[DT$super==se,],GRanges(seqnames=Rle(sub('^chr','',chr)),ranges=IRanges(start=start,end=end),id=enh_id))
	## next for the pchic obtain the SE tissue context specific interactions
	tissues<-unlist(strsplit(ti,','))
	tidx<-which(rowSums(do.call("cbind",lapply(tissues,function(t){
		pchic[[t]]>pchic.cutoff
	})))>0)
	pchic.gr<-with(pchic[tidx,],GRanges(seqnames=Rle(oeChr),ranges=IRanges(start=oeStart,end=oeEnd),id=oeID,ensg=ensg,name=name,pchic=ti,se.tissue=gsub('\\.csv.','',f)))
	mergeByOverlaps(se.gr,pchic.gr)
}

merge.results.se<-lapply(1:nrow(hnisz),function(i){
	message(hnisz[i,]$pchic.tissue)
	mergeSEandPCHiC(hnisz[i,]$f,hnisz[i,]$pchic.tissue,pchic)
})
save(merge.results.se,file=file.path(data.dir,'hnisz_SE_annotation.RData'))
merge.results.non.se<-lapply(1:nrow(hnisz),function(i){
	message(hnisz[i,]$pchic.tissue)
	mergeSEandPCHiC(hnisz[i,]$f,hnisz[i,]$pchic.tissue,pchic,se=0)
})
save(merge.results.non.se,file=file.path(data.dir,'hnisz_Non_SE_annotation.RData'))
hnisz$se.count<-sapply(merge.results.se,nrow)
hnisz$non.se.count<-sapply(merge.results.non.se,nrow)
write.table(hnisz,file=file.path(data.dir,'hnisz_results.csv'),sep=',',row.names=FALSE,quote=FALSE)
