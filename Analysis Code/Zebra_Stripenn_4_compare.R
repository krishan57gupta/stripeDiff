library(export)
library(Hmisc)
library(ggplot2)
library(ggpubr)

stripeNNDF=data.frame()

size_list=c(100000,50000,25000,10000)[c(1:4)]
type_list=c("IMR90","HUVEC")
id_list=c(paste("chr",1:22,sep=""),"chrX")
folder="~/Revision/test_updated"
home="~/Revision/new"
tt="VC1"
for(binSize in size_list){
  for(id in id_list){
    for(type in type_list){
      if(binSize<1000000){
        res=paste(binSize/1000,"kb",sep="")
      }else{
        res=paste(binSize/1000000,"mb",sep="")
      }
      print(paste(type,id,res,"VC",sep="_"))
      
      fileN=paste(home,"/output_dir_",type,"_",binSize,"_",substr(id,4,nchar(id)),"/result_unfiltered.tsv",sep="")
      if(file.exists(fileN)){
        stripe=read.csv(fileN,sep="\t")
        stripeNNDF=rbind(stripeNNDF,cbind(stripe,"dataType"=rep(type,dim(stripe)[1]),"binSize"=rep(binSize,dim(stripe)[1]),"resolution"=rep(res,dim(stripe)[1])))
      }
    } 
  }
}
write.table(stripeNNDF,paste(home,"/stripeNNDF_",tt,".txt",sep=""))
stripeNNDF=read.table(paste(home,"/stripeNNDF_",tt,".txt",sep=""))
stripeNNDF$chr=paste("chr",stripeNNDF$chr,sep="")
stripeDiffDF=read.table(paste(home,"/../1_13/stripeDiffDF_",tt,".txt",sep=""))
stripeDiffDF=stripeDiffDF[stripeDiffDF$binSize %in% c(100000,50000,25000,10000),]
zebraDF=read.table(paste(home,"/../1_13/ZebraDF_",tt,".txt",sep=""))
zebraDF=zebraDF[!duplicated.data.frame(zebraDF[,c("start1","end1")]),]


femo=read.csv(paste(home,"/../1_13/fimo.tsv",sep=""),sep="\t")
table(femo$sequence_name)
femo=femo[femo$sequence_name %in% id_list,]
table(femo$sequence_name)

femo=femo[femo$q.value<=0.05,]

dff=data.frame()
dff=rbind( cbind("chromosome"=femo$sequence_name,"loc"=as.integer((femo$start+femo$stop)/2),"res"=rep("10kb",dim(femo)[1]),
                 "bin"=rep("10000",dim(femo)[1]),"dataType"=rep("IMR90",dim(femo)[1]),"methods"=rep("fimo",dim(femo)[1])),
           cbind("chromosome"=femo$sequence_name,"loc"=as.integer((femo$start+femo$stop)/2),"res"=rep("25kb",dim(femo)[1]),
                 "bin"=rep("25000",dim(femo)[1]),"dataType"=rep("IMR90",dim(femo)[1]),"methods"=rep("fimo",dim(femo)[1])),
           cbind("chromosome"=femo$sequence_name,"loc"=as.integer((femo$start+femo$stop)/2),"res"=rep("50kb",dim(femo)[1]),
                 "bin"=rep("50000",dim(femo)[1]),"dataType"=rep("IMR90",dim(femo)[1]),"methods"=rep("fimo",dim(femo)[1])),
           cbind("chromosome"=femo$sequence_name,"loc"=as.integer((femo$start+femo$stop)/2),"res"=rep("100kb",dim(femo)[1]),
                 "bin"=rep("100000",dim(femo)[1]),"dataType"=rep("IMR90",dim(femo)[1]),"methods"=rep("fimo",dim(femo)[1])),
           cbind("chromosome"=femo$sequence_name,"loc"=as.integer((femo$start+femo$stop)/2),"res"=rep("10kb",dim(femo)[1]),
                 "bin"=rep("10000",dim(femo)[1]),"dataType"=rep("HUVEC",dim(femo)[1]),"methods"=rep("fimo",dim(femo)[1])),
           cbind("chromosome"=femo$sequence_name,"loc"=as.integer((femo$start+femo$stop)/2),"res"=rep("25kb",dim(femo)[1]),
                 "bin"=rep("25000",dim(femo)[1]),"dataType"=rep("HUVEC",dim(femo)[1]),"methods"=rep("fimo",dim(femo)[1])),
           cbind("chromosome"=femo$sequence_name,"loc"=as.integer((femo$start+femo$stop)/2),"res"=rep("50kb",dim(femo)[1]),
                 "bin"=rep("50000",dim(femo)[1]),"dataType"=rep("HUVEC",dim(femo)[1]),"methods"=rep("fimo",dim(femo)[1])),
           cbind("chromosome"=femo$sequence_name,"loc"=as.integer((femo$start+femo$stop)/2),"res"=rep("100kb",dim(femo)[1]),
                 "bin"=rep("100000",dim(femo)[1]),"dataType"=rep("HUVEC",dim(femo)[1]),"methods"=rep("fimo",dim(femo)[1])),
           cbind("chromosome"=zebraDF$chr1,"loc"=as.integer((zebraDF$start1+zebraDF$end1)/2),"res"=zebraDF$resolution,
                 "bin"=zebraDF$binSize,"dataType"=zebraDF$dataType,"methods"=rep("zebra",dim(zebraDF)[1])),
           cbind("chromosome"=stripeNNDF$chr,"loc"=as.integer((stripeNNDF$pos1+stripeNNDF$pos2)/2),"res"=stripeNNDF$resolution,
                 "bin"=stripeNNDF$binSize,"dataType"=stripeNNDF$dataType,"methods"=rep("stripeNN",dim(stripeNNDF)[1])),
           cbind("chromosome"=stripeDiffDF$chrom,"loc"=as.integer((as.integer(stripeDiffDF$upPeak.loc)+as.integer(stripeDiffDF$downPeak.loc))/2*stripeDiffDF$binSize),
                 "res"=stripeDiffDF$resolution,"bin"=stripeDiffDF$binSize,"dataType"=stripeDiffDF$dataType,"methods"=rep("stripeDiff",dim(stripeDiffDF)[1])) )
dff=as.matrix(dff)
dff=data.frame(dff)
sum(is.na(dff))
dff=na.omit(dff)
sum(is.na(dff))
dff$loc=as.numeric(dff$loc)
#########################################################################################################
dfFIMO=dff[dff$methods=="fimo",]
dff12=dff[dff$methods=="stripeNN" | dff$methods=="zebra" | dff$methods=="stripeDiff",]
dff12["nearestCTCFmotif"]=rep(0,dim(dff12)[1])
dff12["nearestCTCFmotifCount"]=rep(0,dim(dff12)[1])

for(loc in 1:dim(dff12)[1]){
  print(paste(dim(dff12)[1],loc,sep="_"))
  dff12$nearestCTCFmotif[loc]=unique(min(abs(dfFIMO$loc[dfFIMO$chromosome==dff12$chromosome[loc]]-dff12$loc[loc])))
  dff12$nearestCTCFmotifCount[loc]=dff12$nearestCTCFmotifCount[loc]+
    sum(abs(dfFIMO$loc[dfFIMO$chromosome==dff12$chromosome[loc]]-dff12$loc[loc])<(as.numeric(dff12$bin[loc])/2))
}

write.table(dff12,paste(home,"/dff12_",tt,".txt",sep=""))
dff12=read.table(paste(home,"/dff12_",tt,".txt",sep=""))

dff12=dff12[dff12$bin %in% c(50000,25000,10000),]

pval=c()
for(i in c("50kb","25kb","10kb")){
  for(j in c("HUVEC","IMR90")){
    for(k in c("stripeNN","zebra")){
      print(paste(i,j,k,sep="_"))
      pval=rbind(pval,cbind(i,j,k,wilcox.test(dff12[dff12$res==i & dff12$dataType==j & dff12$methods=='stripeDiff','nearestCTCFmotifCount'],
                                   dff12[dff12$res==i & dff12$dataType==j & dff12$methods==k,'nearestCTCFmotifCount'],
                                   alternative ='greater' )$p.value))
      
    }
  }
}
colnames(pval)<-c("res","data","method","pval")
pval

dff13=c()
for(i in unique(dff12$dataType)){
    for(k in unique(dff12$res)){
      for(l in unique(dff12$methods)){
        temp=dff12[dff12$dataType==i & dff12$res==k & dff12$methods==l,"nearestCTCFmotifCount"]
        temp1=round(mean(temp),2)
        temp2=round(sd(temp)/sqrt(length(temp)),2)
        dff13=rbind(dff13,c(i,k,l,temp1,temp2))
        
      }
    }
  }

dff13=data.frame(dff13)
colnames(dff13)=c("dataType","res","methods","motifMean","se")
dff13$motifMean=as.numeric(dff13$motifMean)
dff13$se=as.numeric(dff13$se)
dff13$methods=factor(dff13$methods,c("zebra","stripeNN","stripeDiff"))


# pdf(paste(home,"/stripeDiff_stripeZebraNN_boxPlot_",tt,"_3.pdf",sep=""),width=10,height = 10)
plot=ggplot(dff13, aes(x=res, y=motifMean, fill=methods)) +
  geom_bar(stat="identity", position=position_dodge(1))+
  geom_errorbar( aes(x=res, ymin=motifMean-se, ymax=motifMean+se),stat="identity", width=0.5, size=.8,position=position_dodge(1))+
  theme_classic()+ ggtitle(paste("",sep=""))+
  scale_x_discrete(name="resolution", breaks=waiver(), labels=waiver(), limits=c("10kb","25kb","50kb","100kb")[1:3])+
  scale_y_continuous(name="CTCF binding motifs count per stripe", breaks=waiver(), labels=waiver(), limits=NULL)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size=10,family = "Helvetica"))+
  theme(legend.position="top",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))# +coord_flip()
plot+facet_grid(~dataType)
ggsave(filename=paste(home,"/stripeDiff_stripeZebraNN_boxPlot_",tt,"_3.pdf",sep=""),dpi = 300,width=5,height=5)

dev.off()
  