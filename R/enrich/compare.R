library(snpStats)
library(data.table)


AI_DIR<-'/home/ob219/scratch/pid/all_ai_exons_by_chr/snpMatrix_pid_only/'
HT_DIR<-'/home/ob219/scratch/pid/all_ht_exons_by_chr/snpMatrix_pid_only/'

getCADD<-function(file){
	message(file)
	obj<-get(load(file))
	s<-col.summary(obj$gt)
	## get a list of SNPs found in the pid cohort
	#snp.idx<-which(s$RAF > 0 & s$RAF<0.001)
	snp.idx<-which(s$RAF>0)
	i<-obj$info[snp.idx,]
	#as.numeric(i[i$EXON!='',]$CADD_PHRED)
	tmp<-i[grep("splice|frame|start|stop|missense",i$Consequence),.(CADD_PHRED,Gene,row)]
	tmp$CADD_PHRED<-as.numeric(tmp$CADD_PHRED)
	#remove indels
	tmp[tmp$CADD_PHRED!=0,]
}




ai.files<-list.files(path=AI_DIR,pattern='*.RData',full.names=TRUE)
ht.files<-list.files(path=HT_DIR,pattern='*.RData',full.names=TRUE)

ai.cadd<-rbindlist(lapply(ai.files,getCADD))
ht.cadd<-rbindlist(lapply(ht.files,getCADD))

## due to the overlapping nature of genes we have some duplicate variants eliminate these by reselecting by gene list

g<-fread("/home/ob219/scratch/pid/Javierre_prioritised_genes.csv")
g.ai<-unique(subset(g,Trait %in% c('CD','CEL','MS','PBC','RA','SLE','T1D','UC'))$EnsemblID)
g.ht<-unique((subset(g,Trait=='HT')$EnsemblID))

ai.cadd<-ai.cadd[ai.cadd$Gene %in% g.ai]
ht.cadd<-ht.cadd[ht.cadd$Gene %in% g.ht]

ks<-ks.test(ai.cadd$CADD_PHRED,ht.cadd$CADD_PHRED)

performKS<-function(g){
	ks.test(ai.cadd$CADD_PHRED,ht.cadd[ht.cadd$Gene!=g,]$CADD_PHRED)$statistic
}

ksdiff<-do.call('c',lapply(unique(ht.cadd$Gene),performKS))
diff<-abs(ksdiff-ks$statistic)
diff<-data.table(ensg=unique(ht.cadd$Gene),diff=diff)


## let's look at gene size distribution as this could be a confounder larger genes have more chance to have a variant

load("/home/ob219/scratch/pid/gencode.v26.exons.gr.RData")

ht.gr<-subset(e.gr, ensg %in% unique(ht.cadd$Gene))
ht.widths<-sapply(split(ht.gr,ht.gr$ensg),function(gr){
	sum(width(reduce(gr)))
})

## INVESTIGATE WHY SOME OF THESE HAVE HAVE ZERO WIDTH SHOULD NOT BE POSSIBLE
idx<-match(names(widths),diff$ensg)
diff[idx,]$length<-widths
ggplot(diff,aes(x=length,y=log10(diff))) + geom_point()

## looks as if there is a difference but that the difference is the opposite way to that which we suspect
# this could be a function of my gene lists or 
f<-data.table(rbind(data.table(trait='AI',CADD=ai.cadd$CADD_PHRED),data.table(trait='HT',CADD=ht.cadd$CADD_PHRED)))

library(ggplot2)
library(plyr)

caddy <- ddply(f, .(trait), summarize,
                            CADD = unique(CADD),
                            ecdf = ecdf(CADD)(unique(CADD)))

ggplot(caddy, aes(CADD, ecdf, color = trait)) + geom_step()
