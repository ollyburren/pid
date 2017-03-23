library(data.table)
library(GenomicRanges)
library(snpStats)

tabix_bin<-'/home/ob219/bin/htslib/tabix'

## format of this file is that each variant has much data associated with it.
## these routines allow the deconstruction of the annotation data

processVEP<-function(v,vep.header){
	ret<-do.call('rbind',lapply(1:length(v),function(x){
		tmp<-strsplit(v[x],'|',fixed=TRUE)[[1]]
		m<-matrix(tmp,ncol=length(vep.header)-1,byrow=TRUE)
		cbind(m,x)
	}))
	ret<-data.table(ret)
	setnames(ret,c(head(vep.header,-1),'row'))
	ret
}

createInfoDT<-function(str){
	blah<-lapply(strsplit(str,";",fixed=TRUE),function(x) {
		tmp<-strsplit(x,"=",fixed=TRUE)
		ret<-lapply(tmp,'[[',2)
		names(ret)<-sapply(tmp,'[[',1)
		#ret<-data.table(data.frame(ret))
		ret
	})
	## get a list of possible headers
	pos.header<-unique(do.call('c',sapply(blah,names)))
	tlist<-list()
	for(n in pos.header){
		tmp<-sapply(blah,'[[',n)
		if(class(tmp)=='list'){
			tmp<-sapply(tmp,function(x){
				if(length(x)==0){
					return(NA)
				}
				return(x)
			})
			## annoyingly the odd variant is not annotated
			if(n != 'ANN')
				tmp<-as.numeric(tmp)
		}
		tlist[[n]]<-tmp
	}
	#tlist[['ANN']]<-NULL
	tlist<-data.table(data.frame(tlist,stringsAsFactors=FALSE))
	tlist
}


region<-'17:76126851-76139049'
individuals<-c('F000953','F008870')
args<-list(egrep="^#|HGMD_CLASS=DM[^?]|splice|frame|start|stop|missense")

getPossDamSNPs<-function(region,individuals,...){
	chr<-sub("([^:]+):.*","\\1",region)
	start<-as.numeric(sub("[^:]+:([^\\-]+)\\-.*","\\1",region))
	end<-as.numeric(sub("[^:]+:[^\\-]+\\-(.*)","\\1",region))
	#args<-list(...)	
	## load in the merged vcf file.
	vcf.dir<-'/scratch/WGS10K/data/release/latest/merged-vcf/no_hgmd/gnomad/'
	data.dir<-'/scratch/ob219/pid/'
	## load in the file of bridge pid cases
	id.file<-fread(file.path(data.dir,'BRIDGE-PID_20161012release_846_index_cases.csv'))
	setnames(id.file,make.names(names(id.file)))
	#fname<-'chr%s_agg3_dedup_vep.vcf.gz'
	fname<-'chr%s_agg3_dedup_vep.vcf.gz'
	vf<-file.path(vcf.dir,sprintf(fname,chr))
	my.pipe<-pipe(paste(tabix_bin,'-H',vf,region))
	header<-tail(scan(my.pipe,what=character(),sep="\n",quiet=TRUE),n=1)
	close(my.pipe)
	cnames<-unlist(strsplit(header,"\t"))
	## variant annotation is another thing to capture
	my.pipe<-pipe(paste(tabix_bin,'-H',vf,"| egrep  '^##INFO=<ID=ANN' | sed -e s'/^.*Format: \\([^\\][^\\]*\\)\\.*/\\1/'"))
	vep<-scan(my.pipe,what=character(),sep="\n",quiet=TRUE)
	vep<-gsub('"','',vep)
	vep.header<-strsplit(vep,'\\|')[[1]]
	cols.want<-paste(c(1:10,which(cnames %in% id.file$WGS.ID)),sep=",",collapse=",")
	cmd_mod<-sprintf("| cut -f%s",cols.want)
	if("egrep" %in% names(args)){
		egrep_cmd<-sprintf('| egrep "%s"',args[['egrep']])
		cmd_mod<-paste(egrep_cmd,cmd_mod)
	}
	## this is slow reading vcf into R
	overall_cmd<-paste(tabix_bin,vf,region,cmd_mod)
	message(overall_cmd)
	tmp<-as.data.frame(fread(overall_cmd,sep="\t",header=FALSE,stringsAsFactors=FALSE))
	colnames(tmp)<-cnames[c(1:10,which(cnames %in% id.file$WGS.ID))]
	gt<-tmp[,10:ncol(tmp)]
	if(nrow(gt)==0)
	  return(NA)
	info<-tmp[,1:9]
	sm<-apply(gt,1,function(x) sub("0\\/0.*","1",x))
	sm<-apply(sm,1,function(x) sub("(0\\/1).*|(1\\/0).*","2",x))
	sm<-apply(sm,1,function(x) sub("1\\/1.*","3",x))
	sm<-t(apply(sm,1,function(x) as.raw(sub("[0-9]\\/[0-9]","0",x))))
	colnames(sm)<-1:nrow(info)
	rownames(sm)<-colnames(gt)
	sm<-new("SnpMatrix", sm)
	info.dt<-createInfoDT(info$INFO)
	obj<-list(map=info,gt=sm,info=info.dt)
	#get VEP stuff - note that we get multiple entries per SNP.
	vep<-processVEP(obj$info$ANN,vep.header)
	## get SNPs that have a high chance to mess up the protein
	interesting.vep<-vep[grep(args$egrep,vep$Consequence),c('Consequence','CADD_PHRED','SIFT','aaalt','aaref','PolyPhen','EXON','ExAC_AF','row'),with=FALSE]
	#c.idx<-as.numeric(unique(subset(vep,IMPACT %in% c('HIGH','MODERATE') & CANONICAL=='YES' & BIOTYPE=='protein_coding')$row))
	## use snpMatrix to find out which ones have variation
	#csum<-col.summary(obj$gt[,c.idx])
	csum<-col.summary(obj$gt[individuals,as.numeric(unique(interesting.vep$row))])
	idx<-as.numeric(rownames(csum[which(csum$RAF>0),]))
	gt.res<-gt[idx,individuals]
	if(length(individuals)==1){
		## in this case we get a scalar 
		gt.res<-data.frame(idx,gt.res)
		names(gt.res)<-c('id',individuals)
	}else{
		gt.res$id<-rownames(gt.res)
	}
	merge(gt.res,interesting.vep,by.x='id',by.y='row')
}

## get a region that we are interested in 
## get bridge identifiers that we are intested in
#region<-'1:114356433-114414381'
## load in the list of genes we wish to check
(load("/scratch/ob219/pid/loss_merge_genes/pid.se.RData"))
## need to use biomart to get gene locations 
library(biomaRt)

mart <- useMart('ENSEMBL_MART_ENSEMBL',host="grch37.ensembl.org")
ensembl.gene <- useDataset("hsapiens_gene_ensembl",mart=mart)

gene.details<-getBM(
                filters= c("ensembl_gene_id"),
                attributes= c('ensembl_gene_id','chromosome_name','start_position','end_position'),
                values= list(ensembl_gene_id=out.tab$p.ensg),
                mart= ensembl.gene)
DT<-with(gene.details,data.table(ensg=ensembl_gene_id,region=sprintf("%s:%s-%s",chromosome_name,start_position,end_position)))
setkey(DT,ensg)
setkey(out.tab,p.ensg)
ot<-out.tab[DT]
by.reg<-split(ot$sample,ot$region)

#for(reg in names(by.reg)){
#	message(sprintf("Doing %s",reg))
#	ind<-by.reg[[reg]]
#	getPossDamSNPs(reg,ind)
#}
out<-lapply(seq_along(by.reg),function(x){
	reg<-names(by.reg)[x]
	message(sprintf("Doing %s",reg))
        ind<-by.reg[[reg]]
        getPossDamSNPs(reg,ind,egrep="^#|HGMD_CLASS=DM[^?]|splice|frame|start|stop|missense")
})

names(out)<-names(by.reg)
## TMC6/8 potentially interesting 
## TCF6
## PLCG2
## NOD2
## look in more detail next week
save(out,file="/scratch/ob219/pid/loss_merge_genes/pid_se_possible_lof.RData")
