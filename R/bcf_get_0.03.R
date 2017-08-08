library(data.table)
library(GenomicRanges)
library(snpStats)

bcftools_bin<-'~/bin/bcftools-1.4/bcftools'

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
	pos.header<-unique(c(do.call('c',lapply(blah,names)),'GNOMAD_AF'))
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

robustDTfread<-function(cmd){
	tryCatch(as.data.frame(fread(cmd,sep="\t",header=FALSE,stringsAsFactors=FALSE)),warning=function(w) {print(sprintf("Warning=%s fread CMD=%s",w,cmd))},error=function(e){print(sprintf("Error=%s fread CMD=%s",e,cmd));return(NA)})
}

robustDTfread<-function(cmd){
	tryCatch(as.data.frame(fread(cmd,sep="\t",header=FALSE,stringsAsFactors=FALSE)),error=function(e){print(sprintf("Error=%s fread CMD=%s",e,cmd));return(NA)})
}
## this is prototype command for bcftools
## finish converting to BCFTools
## look at genes outside of the PID list
## reread the bioarxiv paper on reg variation
## reread BeviMed paper.
getPossDamSNPs<-function(region,individuals,...){
	ar<-list(...)
	chr<-sub("([^:]+):.*","\\1",region)
	start<-as.numeric(sub("[^:]+:([^\\-]+)\\-.*","\\1",region))
	end<-as.numeric(sub("[^:]+:[^\\-]+\\-(.*)","\\1",region))
	ind<-paste(individuals,sep=',',collapse=',')
	#args<-list(...)
	## load in the merged vcf file.
	vcf.dir<-'/scratch/WGS10K/data/release/latest/merged-vcf/no_hgmd/gnomad/'
	data.dir<-'/scratch/ob219/pid/'
	fname<-'chr%s_agg3_dedup_vep_gnomad.bcf'
	vf<-file.path(vcf.dir,sprintf(fname,chr))
	header_cmd<-sprintf("%s view -s %s --header-only %s",bcftools_bin,ind,vf)
	my.pipe<-pipe(header_cmd)
	header<-tail(scan(my.pipe,what=character(),sep="\n",quiet=TRUE),n=1)
	close(my.pipe)
	cnames<-unlist(strsplit(header,"\t"))
	## should be able to do in one pass but this is so quick !
	my.pipe<-pipe(sprintf("%s | egrep -e '%s' | sed -e s'%s'",header_cmd,'^##INFO=<ID=ANN','/^.*Format: \\([^\\][^\\]*\\)\\.*/\\1/'))
	## variant annotation header is another thing to capture
	vep<-scan(my.pipe,what=character(),sep="\n",quiet=TRUE)
	close(my.pipe)
	vep<-gsub('"','',vep)
	vep.header<-strsplit(vep,'\\|')[[1]]
	if("egrep" %in% names(ar)){
		body_cmd<-sprintf("%s view %s -r %s -s %s -Ov -H | egrep -e '%s'",bcftools_bin,vf,region,ind,ar[['egrep']])
	}else{
		body_cmd<-sprintf("%s view %s -r %s -s %s -Ov -H",bcftools_bin,vf,region,ind)
	}
	message(body_cmd)
	#tmp<-as.data.frame(fread(body_cmd,sep="\t",header=FALSE,stringsAsFactors=FALSE))
	tmp<-robustDTfread(body_cmd)
	if(is.na(tmp))
		return(NA)
	colnames(tmp)<-cnames
	gt<-tmp[,10:ncol(tmp)]
	if(length(individuals)==1){
		## add dummy individual which we can remove at the end
		gt<-cbind(gt,DUMMY=rep('0/0',length(gt)))
		colnames(gt)<-c(individuals,'DUMMY')
	}
	if(nrow(gt)==0)
	  return(NA)
	info<-tmp[,1:9]
	sm<-apply(gt,1,function(x) sub("0\\/0.*","1",x))
	sm<-as.matrix(apply(sm,1,function(x) sub("(0\\/1).*|(1\\/0).*","2",x)))
	sm<-as.matrix(apply(sm,1,function(x) sub("1\\/1.*","3",x)))
	sm<-t(apply(sm,1,function(x) as.raw(sub("[0-9]\\/[0-9]","0",x))))
	## single snps don't need to be transposed
	message(length(colnames(sm)))
	if(length(colnames(sm))!=0)
		sm<-t(sm)
	colnames(sm)<-1:nrow(info)
	rownames(sm)<-colnames(gt)
	sm<-new("SnpMatrix", sm)
	info.dt<-createInfoDT(info$INFO)




	## possibly only interested in GNOMAD_AF to start with
	obj<-list(map=info,gt=sm,info=info.dt)
	#get VEP stuff - note that we get multiple entries per SNP.
	vep<-processVEP(obj$info$ANN,vep.header)
	## get SNPs that have a high chance to mess up the protein
	interesting.vep<-vep[grep(ar$egrep,vep$Consequence),c('SYMBOL','Gene','BIOTYPE','Existing_variation','Consequence','CADD_PHRED','SIFT','Amino_acids','PolyPhen','EXON','row'),with=FALSE]
	## get info on GNOMAD allele freq
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
	## GNOMAD AF
	gaf<-data.frame(id=idx,GNOMAD_AF=info.dt[idx,]$GNOMAD_AF)
	tmp<-merge(gt.res,interesting.vep,by.x='id',by.y='row')
	merge(tmp,gaf,by.x='id',by.y='id')
}


(load("/scratch/ob219/pid/loss_merge_genes/all.se.0.03.RData"))
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
out<-lapply(seq_along(by.reg),function(x){
	reg<-names(by.reg)[x]
	message(sprintf("Doing %s",reg))
        ind<-unique(by.reg[[reg]])
        tmp<-getPossDamSNPs(reg,ind,egrep="^#|HGMD_CLASS=DM[^?]|splice|frame|start|stop|missense")
})

names(out)<-names(by.reg)

## TMC6/8 potentially interesting
## TCF6
## PLCG2
## NOD2
## look in more detail next week
save(out,file="/scratch/ob219/pid/loss_merge_genes/all_se_possible_lof0.03.RData")
stop();
load("/scratch/ob219/pid/loss_merge_genes/all_se_possible_lof.RData")

## remove those for which there are no variants

out.f<-out[sapply(out,is.data.frame)]
## only look at protein_coding genes for time being
bl<-sapply(out.f,function(x){
	if(any(x$BIOTYPE=='protein_coding'))
		return(TRUE)
	return(FALSE)
})
out.f<-out.f[bl]
## next get those with CADD scores above 15
out.f<-sapply(out.f,function(x){
	x$CADD_PHRED<-as.numeric(x$CADD_PHRED)
	x$GNOMAD_AF<-as.numeric(x$GNOMAD_AF)
	x[!is.na(x$CADD_PHRED) & x$CADD_PHRED>15 & (x$GNOMAD_AF<0.05 | is.na(x$GNOMAD_AF)),]
})
out.f<-out.f[sapply(out.f,nrow)!=0]
unique(do.call('c',lapply(out.f,function(x) unique(x$SYMBOL))))
rbindlist(out.f)
## ok fi
