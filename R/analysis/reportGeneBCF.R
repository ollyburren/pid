library(data.table)



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

## Special for just looking at PIK3C2B

parseVCF<-function(){
	vf<-'/home/ob219/scratch/pid/VCF/tmp/PIK3C2B.bcf'
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
	body_cmd<-sprintf("%s view %s -Ov -H",bcftools_bin,vf)
	#message(body_cmd)
	tmp<-robustDTfread(body_cmd)
	if(any(is.na(tmp)))
		return(NA)
	colnames(tmp)<-cnames
	info<-tmp[,1:9]
	info.dt<-createInfoDT(info$INFO)
	vep<-processVEP(info.dt$ANN,vep.header)
	all.info<-merge(info.dt,vep,by.x='row',by.y='row')[,.(WGS10K_AF,GNOMAD_AF,OPR,SYMBOL,Gene,BIOTYPE,Existing_variation,Consequence,CADD_PHRED,Amino_acids,row)]
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
parseVCF()
