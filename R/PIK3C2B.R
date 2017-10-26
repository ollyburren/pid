library(data.table)
library(biomaRt)


processVEP<-function(v,vep.header){
	ret<-do.call('rbind',lapply(1:length(v),function(x){
		tmp<-strsplit(v[x],'|',fixed=TRUE)[[1]]
		m<-matrix(tmp,ncol=length(vep.header)-1,byrow=TRUE)
		cbind(m,x)
	}))
	ret<-data.table(ret)
	setnames(ret,c(head(vep.header,-1),'row'))
	ret$row<-as.numeric(ret$row)
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
	tlist$row<-1:nrow(tlist)
	tlist
}

robustDTfread<-function(cmd){
	tryCatch(as.data.frame(fread(cmd,sep="\t",header=FALSE,stringsAsFactors=FALSE)),warning=function(w) {print(sprintf("Warning=%s fread CMD=%s",w,cmd))},error=function(e){print(sprintf("Error=%s fread CMD=%s",e,cmd));return(NA)})
}

robustDTfread<-function(cmd){
	tryCatch(as.data.frame(fread(cmd,sep="\t",header=FALSE,stringsAsFactors=FALSE)),error=function(e){print(sprintf("Error=%s fread CMD=%s",e,cmd));return(NA)})
}



# vf - vcf file
# bf - bed file of specific regions to extract

parseVCF<-function(vf,bf){
  message(sprintf("Processing %s",bf))
	if(!file.exists(vf))
		return(NA)
	header_cmd<-sprintf("%s view --header-only %s",bcftools_bin,vf)
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
	body_cmd<-sprintf("%s view %s -Ov -H -R %s",bcftools_bin,vf,bf)
	#message(body_cmd)
	tmp<-robustDTfread(body_cmd)
	if(any(is.na(tmp)))
		return(NA)
	colnames(tmp)<-cnames
	info<-tmp[,1:9]
	info.dt<-cbind(createInfoDT(info$INFO),info[,c('POS','REF','ALT')])
	vep<-processVEP(info.dt$ANN,vep.header)
	all.info<-merge(info.dt,vep,by.x='row',by.y='row')[,.(POS,REF,ALT,WGS10K_AF,GNOMAD_AF,OPR,SYMBOL,Gene,BIOTYPE,Existing_variation,Consequence,CADD_PHRED,Amino_acids,row)]
	## next remove those that cannot be loss of function
	#all.info<-all.info[grep("HGMD_CLASS=DM[^?]|splice|frame|start|stop|missense",all.info$Consequence),]
	inds<-names(tmp)[10:length(names(tmp))]
	stub<-rbindlist(lapply(inds,function(i){
			igt<-tmp[[i]]
			igt<-sub("0\\/0.*","1",igt)
			igt<-sub("(0\\/1).*|(1\\/0).*","2",igt)
			igt<-sub("1\\/1.*","3",igt)
			igt<-sub("[0-9]\\/[0-9]","0",igt)
			data.table(ind=i,gt=igt,row=(1:length(igt)))
	}))
	merge(stub,all.info,by.x='row',by.y='row',allow.cartesian=TRUE)
}

## script to extract PIK3C2B coding variants across all cohorts

bcftools_bin<-'~/bin/bcftools-1.4/bcftools'
## manifest file
m_file<-'/scratch/WGS10K/data/projects/BRIDGE_5K_PAPER/data_table_20170822.txt'
mani<-fread(m_file)
setnames(mani,make.names(names(mani)))
gene<-'PIK3C2B'
out.dir<-file.path('/scratch/ob219/pid/ANNOTATIONS/',gene)
if(!file.exists(out.dir))
  dir.create(out.dir)
## get exonic coordinates

ensembl_archive <- 'feb2014.archive.ensembl.org'
ensembl <- useMart(host=ensembl_archive, biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
gene.details<-getBM(c("ensembl_gene_id","gene_biotype","external_gene_id","hgnc_symbol", "chromosome_name","start_position", "end_position"),
  mart=ensembl,filters=list(hgnc_symbol=gene))

## next extract the relevant exon coords from the correct file

efile<-sprintf("/home/ob219/scratch/pid/ANNOTATIONS/exons_by_chr/chr%s_gencode.v26lift37.bed",gene.details[1,]$chromosome_name)
ofile<-file.path(out.dir,'exons.bed')
cmd<-sprintf("grep %s %s > %s",gene,efile,ofile)
system(cmd)

## get the data we want in R

vcf_file <- sprintf("/scratch/WGS10K/data/release/20170104-A/merged-vcf/no_hgmd/gnomad/chr%s_agg3_dedup_vep_gnomad.bcf",gene.details[1,]$chromosome_name)

all.cohort.snps<-parseVCF(vcf_file,ofile)[Gene==gene.details[1,]$ensembl_gene_id,]
## next we need to filter so just have the gene we want

## huge data.table

## add an id based on chr position and alleles

all.cohort.snps[,id:=paste(1,POS,REF,ALT,sep=':')]
setkey(all.cohort.snps,ind)
setkey(mani,X.WGS_ID)

full<-merge(all.cohort.snps,mani,by.x='ind',by.y='X.WGS_ID')

## let's get allele counts

## next merge in manifest file

fui<-unique(full,by=c('ind','id'))

## for counts create a table of unique counts

## perhaps limit to LOF variants

all.cons<-unique(full$Consequence)

lof.var<-c('stop_gained','frameshift_variant','missense_variant',all.cons[grep('splice',all.cons)])

#lof.var<-c('stop_gained','frameshift_variant','missense')

ac<-full[AFFECTED=='Y' & as.numeric(OPR) > 0.98 & Consequence %in% lof.var,list(AC=.N),by=c('id','PROJECT','gt')]

## compute allele frequencies
ac[,a2.f:=(AC * (as.numeric(gt) -1))/7226]
overall.ac <- ac[,list(total=sum(AC)),by=c('id','gt')]
overall.af <- overall.ac[,af:=(total * (as.numeric(gt) -1))/(7226*2)][,list(ol.af=sum(af)),by=id]
overall.af[,maf:=ol.af]
overall.af[maf>0.5,maf:=1-maf]


getUniqueByCohort<-function(cname){
	bycohort<-ac[,list(tot=sum(a2.f)),by=c('id','PROJECT')]
	bycohort[,byc:=cname]
	bycohort[PROJECT!=cname,byc:='other']
	bycohort <- bycohort[,list(T=sum(tot)),by=c('id','byc')]
	overall.ac <- ac[,list(total=sum(AC)),by=c('id','gt')]
	overall.af <- overall.ac[,af:=(total * (as.numeric(gt) -1))/(7226*2)][,list(ol.af=sum(af)),by=id]
	overall.af[,maf:=ol.af]
	overall.af[maf>0.5,maf:=1-maf]
	##  next add some ANNOTATIONS
	# is variant only found in PID ?
	tmp<-overall.af[id %in% bycohort[byc=='other' & T==0,]$id & maf>0,]
	tmp[,cohort:=cname]
}

all.cohorts<-unique(full$PROJECT)

res<-rbindlist(lapply(all.cohorts,getUniqueByCohort))
table(res$cohort)

getUniqueByCohort('PID')

bycohort[,bypid:='pid']
bycohort[PROJECT!='PID',bypid:='other']
bycohort <- bycohort[,list(T=sum(tot)),by=c('id','bypid')]

overall.ac <- ac[,list(total=sum(AC)),by=c('id','gt')]
overall.af <- overall.ac[,af:=(total * (as.numeric(gt) -1))/(7226*2)][,list(ol.af=sum(af)),by=id]
overall.af[,maf:=ol.af]
overall.af[maf>0.5,maf:=1-maf]
##  next add some ANNOTATIONS
# is variant only found in PID ?
unique.to.pid<-overall.af[id %in% bycohort[bypid=='other' & T==0,]$id & maf>0,]
## across all cohorts see what the number of unique variants is - so  as above but replace PID
## with other cohort (create a function to do this)

##

acm<-data.table::melt(ac,id.vars=c('id','PROJECT','gt'),measured.vars='AC')
dcast(acm,id~PROJECT+gt+variable)
