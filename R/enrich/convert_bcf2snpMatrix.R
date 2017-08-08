library(data.table)
library(snpStats)
library(optparse)
bcftools_bin<-'~/bin/bcftools-1.4/bcftools'
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
        tryCatch(as.data.frame(fread(cmd,sep="\t",header=FALSE,stringsAsFactors=FALSE)),error=function(e){print(sprintf("Error=%s fread CMD=%s",e,cmd));return(NA)})
}

saveBCF2SM<-function(vf='/home/ob219/scratch/pid/all_ai_exons_by_chr/bcf/ai_18.bcf',ofile='/home/ob219/scratch/pid/all_ai_exons_by_chr/snpMatrix/ai_18.RData'){
	## for time being don't sweat the individuals
	#vf<-'/home/ob219/scratch/pid/all_ai_exons_by_chr/bcf/ai_18.bcf'
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
	tmp<-robustDTfread(body_cmd)
	colnames(tmp)<-cnames
	info<-tmp[,1:9]
	gt<-tmp[,10:ncol(tmp)]
	sm<-apply(gt,1,function(x) sub("0\\/0.*","1",x))
	sm<-as.matrix(apply(sm,1,function(x) sub("(0\\/1).*|(1\\/0).*","2",x)))
	sm<-as.matrix(apply(sm,1,function(x) sub("1\\/1.*","3",x)))
	sm<-t(apply(sm,1,function(x) as.raw(sub("[0-9]\\/[0-9]","0",x))))
	## single snps don't need to be transposed
	colnames(sm)<-1:nrow(info)
	rownames(sm)<-colnames(gt)
	sm<-new("SnpMatrix", sm)
	info.dt<-createInfoDT(info$INFO)
	#get VEP stuff - note that we get multiple entries per SNP.
	vep<-processVEP(info.dt$ANN,vep.header)
	## get values from VEP with the best CADD score
	## for this to work we need to make "" be zero (this score should be impossible) - I think that these are INDELS for which CADD is not available
	vep[vep$CADD_PHRED=="",]$CADD_PHRED<-"0"
	cadd.max<-vep[vep[,.I[which.max(CADD_PHRED)],by=row]$V1,]
	cadd.max<-cadd.max[,.(Consequence,CADD_PHRED,EXON,CCDS,Existing_variation,Gene,SYMBOL,row)]
	if(nrow(cadd.max) != length(colnames(sm)))
		message("Warning number of SNPs and annotations disagree")
	obj<-list(map=info,gt=sm,info=cadd.max)
	message(sprintf("Saving %d SNPs to %s",nrow(cadd.max),ofile))
	save(obj,file=ofile)
}

option_list = list(
	make_option(c("-f", "--file"), type="character", default=NULL, 
              help="vcf file to convert", metavar="character"),
	make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 
## option parsing 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$f)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
message(sprintf("Converting %s to out file %s",opt$file,opt$out))
## convert
saveBCF2SM(vf=opt$file,ofile=opt$out)


