library(data.table)
library(GenomicRanges)
library(snpStats)

## annotate each CNV event with overlap with various categories


DATA.DIR<-'/home/ob219/scratch/pid/ANNOTATIONS/'

## promoter regions
prom<-readRDS(file.path(DATA.DIR,'promoters.RDS'))

## exon regions
exon<-readRDS(file.path(DATA.DIR,'exons.RDS'))

## pchic superenhancer (for time being)
pse<-readRDS(file.path(DATA.DIR,'pchic_hnisz','se.gr.RDS'))

genes<-readRDS(file.path(DATA.DIR,'genes.RDS'))

## create a table of CNV's and at this stage annotate true or false for whether there is an overlap for one of the above features

## load in all the events (bnd will need to be handled differently)

manta.files<-list.files(path='/home/ob219/scratch/pid/MANTA/',pattern="*.RDS",full.names=TRUE)
mevents<-lapply(manta.files,readRDS)
names(mevents)<-gsub("pid\\_([^_]+).*","\\1",basename(manta.files))
## drop bnd for time being
mevents<-mevents[names(mevents) != 'bnd']

computeFeatOl<-function(e.gr,f.gr){
  out<-logical(length=length(e.gr))
  ol<-unique(as.matrix(findOverlaps(e.gr,f.gr))[,1])
  out[ol]<-TRUE
  out
}

annoMatrix<-function(e){
  tmp<-do.call('cbind',lapply(list(prom=prom,exon=exon,pse=pse),function(f.gr){
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

## ok res contains the 'interesting' results where an event overlaps at least one feature.
## next we want to make a more detailed annotation of what these features are by sample


f2s<-function(feat){
  flist<-list(prom=prom,exon=exon,pse=pse)
  rbindlist(lapply(names(flist),function(n){
    f.gr<-flist[[n]]
    of<-as.matrix(findOverlaps(f.gr,feat))
    of<-split(of[,2],of[,1])
    #make things quicker
    #feat$eid<-paste(feat$bridge_id,feat$uid,feat$gt,sep=':')
    m<-mergeByOverlaps(feat,f.gr)
    data.table(gene=m$ensg,type=n,bridge_id=m$bridge_id,gt=m$gt,eid=m$uid)
    #lapply(split(m$eid,m$ensg),unique)
  }))
}

fr<-rbindlist(lapply(names(res),function(n){
    tmp<-f2s(res[[n]])
    tmp$event=n
    tmp
}))

## first combine over event type
all<-unique(fr[,.(gene,bridge_id)])

tmp<-lapply(split(all$bridge_id,all$gene),unique)

## ok add in the gene coords for each and then the samples
bcftools_bin<-'~/bin/bcftools-1.4/bcftools'
cmd<-"%s view %s -r %s -s %s -Ob | %s view - -i 'AC!=0' -Ob -o %s"
lof.dir<-'/home/ob219/scratch/pid/VCF/exonic_only/no_monomorphs/lof/'
out.dir<-'/home/ob219/scratch/pid/VCF/exonic_only/no_monomorphs/lof/gene_events/'
genes.by.ensg<-split(genes,genes$ensg)
all.commands<-sapply(names(tmp),function(n){
  g<-genes.by.ensg[[n]]
  if(length(g)!=1)
    return(NA)
  region<-sprintf("%s:%d-%d",seqnames(g),start(g),end(g))
  samples<-paste(tmp[[n]],sep=',',collapse=',')
  outfile<-file.path(out.dir,sprintf("%s_lof.bcf",n))
  vcf.file=file.path(lof.dir,sprintf("chr%s_gencode.v26lift37.bcf",seqnames(g)))
  sprintf(cmd,bcftools_bin,vcf.file,region,samples,bcftools_bin,outfile)
})

write(all.commands[!is.na(all.commands)],file='/scratch/ob219/pid/VCF/gene_events_lof.txt')

## run these on the q - the ones that are NA are the non protein coding I think.

## want to parse all the results and create a simplified dossier for each gene that we can refer to
## info to capture is allele freq, and samples and alleles




## what we want is the coords for each gene linked to a list of samples - we can then use bcftools to get a vcf files for each LoF genes
## chr start stop sample1 sample2 sample3 - then we can write a script to generate bcftools jobs to create lots of small jobs to extract relevant lof SNPs and genotypes

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

all.v.genes<-names(all.commands)[!is.na(all.commands)]

g<-all.v.genes[1]

getShortReport<-function(g){
  vcf.dir<-'/home/ob219/scratch/pid/VCF/exonic_only/no_monomorphs/lof/gene_events/'
  gbcf<-sprintf("%s%s_lof.bcf",vcf.dir,g)
  if(!file.exists(gbcf))
    return(NA)
  header_cmd<-sprintf("%s view --header-only %s",bcftools_bin,gbcf)
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
  body_cmd<-sprintf("%s view %s -Ov -H",bcftools_bin,gbcf)
  tmp<-robustDTfread(body_cmd)
  if(any(is.na(tmp)))
    return(NA)
  colnames(tmp)<-sub("^#","",cnames)
  gt<-tmp[,10:ncol(tmp)]
  if(class(gt)=='character'){
    names(gt)<-rep(cnames[10],length(gt))
    sm<-sub("0\\/0.*","1",gt)
    sm<-sub("(0\\/1).*|(1\\/0).*","2",sm)
    sm<-sub("1\\/1.*","3",sm)
    sm<-t(sub("[0-9]\\/[0-9]","0",sm))
  }else{
    if(nrow(gt)==0)
      return(NA)
    sm<-apply(gt,1,function(x) sub("0\\/0.*","1",x))
    sm<-as.matrix(apply(sm,1,function(x) sub("(0\\/1).*|(1\\/0).*","2",x)))
    sm<-as.matrix(apply(sm,1,function(x) sub("1\\/1.*","3",x)))
    sm<-t(apply(sm,1,function(x) sub("[.0-9][\\/:][.0-9].*","0",x)))
    if(length(colnames(sm))==0)
      sm<-t(sm)
  }
  rel.gt<-apply(sm,1,function(x){
      ## remove where homozygous to reference
      tmp<-x[!x %in% c(0,1)]
      paste(names(tmp),tmp,sep=':',collapse=',')
  })
  info<-tmp[,1:9]
  info.dt<-createInfoDT(info$INFO)
  unique(data.table(gene=g,chr=tmp$CHROM,position=tmp$POS,gnomad_af=info.dt$GNOMAD_AF,wgs10k=info.dt$WGS10K_AF,gt=rel.gt))
}


sres<-lapply(all.v.genes,function(ge){
  message(ge)
  getShortReport(ge)
})

all.sres<-rbindlist(sres[!is.na(sres)])

saveRDS(all.sres,file='/home/ob219/scratch/pid/ANNOTATIONS/lof_snps.RDS')
