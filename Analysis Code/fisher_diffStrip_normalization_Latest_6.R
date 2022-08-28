library(ggplot2)
library(export)
library(Hmisc)
library(doParallel)
library(foreach)
home="~/Revision"

normMethod_list=c("NOnorm_400","KRnorm_400","VCnorm_400","VCSQRTnorm_400","SCNnorm_400")[]
# cl=length(normMethod_list)
# registerDoParallel(cl)

bar_pval_df=data.frame(matrix(ncol = 7, nrow = 0))
bar_pval_df_c<-c("normMethod","res","pval","data","observ","expect","fold")

for(normMethod in normMethod_list){
# foreach::foreach(normMethod=normMethod_list,.combine=cbind)%dopar%{
  for (pval in c(.05)){
    folder=paste(home,"/1_2/data_Latest_6",normMethod,"_",pval,sep="")
    folder2=paste(home,"/1_2/data","_",pval,sep="")
    folder1=paste(home,"/1_2/data",normMethod,sep="")
    if(!dir.exists(folder)){
      dir.create(folder)
    }
    if(!dir.exists(folder1)){
      dir.create(folder1)
    }
    if(!dir.exists(folder2)){
      dir.create(folder2)
    }
    # direction left means, horihental stripe in upper right traingle
    # direction right means, horihental stripe in lower left traingle
    id_list=c(paste("chr",1:22,sep=""),"chrX")
    binSizeList=c(100000,50000,25000,10000)[4:4]
    #foreach::foreach(binSize=binSizeListt,.combine=cbind)%dopar%{
    for(binSize in binSizeList) {
      # binSize=1000000
      if(binSize<1000000){
        res=paste(binSize/1000,"kb",sep="")
      }else{
        res=paste(binSize/1000000,"mb",sep="")
      }
      
      
      cm=matrix(0,ncol=3,nrow=3)
      colnames(cm)=c("IMR90_strengthened","IMR90_weakened","insignificant")
      rownames(cm)=c("IMR90_upregulated","IMR90_downregulated","insignificant")
      
      
      type_list=c("IMR90","HUVEC")
      print(paste(res,pval,normMethod,res,sep="_"))
      dataName=paste("A_6_stripe_",type_list[1],"_",type_list[2],"_",pval,sep="")
      if(file.exists(paste(folder,"/",dataName,"_",res,".rds",sep=""))){
        stripes1=readRDS(paste(folder,"/",dataName,"_",res,".rds",sep=""))
      }else{
        stripes1=list()
      }
      dataName1=paste("gene_",type_list[1],"_",type_list[2],"_",pval,sep="")
      if(file.exists(paste(folder2,"/",dataName1,".rds",sep=""))){
        gene1=readRDS(paste(folder2,"/",dataName1,".rds",sep=""))
        gene1L=dim(gene1)[1]
      }else{
        gene1=c()
        gene1L=1
      }
      dataName2=paste("A_6_FinalResult_",type_list[1],"_",type_list[2],"_",pval,sep="")
      if(file.exists(paste(folder,"/",dataName2,"_",res,".rds",sep=""))){
        finalResult=readRDS(paste(folder,"/",dataName2,"_",res,".rds",sep=""))
      }else{
        finalResult=list()
      }
      
      stripes=c()
      flag=0
      for(id in id_list){
        print(paste("stripe reading ",id,sep=""))
        files=paste(home,"/test_updated/output_updated_",normMethod,"/",type_list[1],"_",type_list[2],"_",
                    id,"_",res,"_obs",normMethod,"/deduplicated_in_",type_list[1],"_not_",type_list[2],"_",id,"_stripes.txt",sep="")
        if(file.exists(files)){
          temp=read.table(files)
          if(flag==0){
            colname=c("in","out","res","comp",temp[1,])
            flag=1
          }
          stripes=rbind(stripes,cbind("in"=rep(type_list[1],dim(temp)[1]-1),"out"=rep(type_list[2],dim(temp)[1]-1),
                                      "res"=rep(res,dim(temp)[1]-1),"comp"=rep(paste(type_list[1],"_",type_list[2],sep=""),dim(temp)[1]-1),temp[-1,]))
        }
        files=paste(home,"/test_updated/output_updated_",normMethod,"/",type_list[1],"_",type_list[2],"_",
                    id,"_",res,"_obs",normMethod,"/deduplicated_in_",type_list[2],"_not_",type_list[1],"_",id,"_stripes.txt",sep="")
        if(file.exists(files)){
          temp=read.table(files)
          stripes=rbind(stripes,cbind("in"=rep(type_list[2],dim(temp)[1]-1),"out"=rep(type_list[1],dim(temp)[1]-1),
                                      "res"=rep(res,dim(temp)[1]-1),"comp"=rep(paste(type_list[1],"_",type_list[2],sep=""),dim(temp)[1]-1),temp[-1,]))
        }
      }
      colnames(stripes)<-colname
      stripes<-na.omit(stripes)
      dim(stripes)
      
      stripes1=readRDS(paste(folder,"/",dataName,"_",res,".rds",sep=""))
      stripes1=as.data.frame(stripes1)
      table(unlist(stripes1$sample1StripeType))
      dim(stripes1)
      
      ################################################################################################################
      gene=read.csv(paste(home,"/1_24/gene_exp.diff",sep=""), sep = "\t")
      # gene=gene[gene$significant=="yes" & gene$status=="OK",]
      # gene=gene[gene$status=="OK",]
      if(gene1L<dim(gene)[1]){
        for(kk in (gene1L+1):dim(gene)[1]){
          print(kk)
          geneLoc=strsplit(gene$locus[kk],":")
          chrL=geneLoc[[1]][1]
          FL=as.numeric(strsplit(geneLoc[[1]][2],"-")[[1]][1])
          SL=as.numeric(strsplit(geneLoc[[1]][2],"-")[[1]][2])
          if(gene[kk,]$p_value<=pval){
            if(gene[kk,]$log2.fold_change>=0){# beacuse first sample is HUVEC,but testing from gene file log2(value_2/value_1) has taken
              gene1=rbind(gene1,cbind(gene[kk,],"regulation"="IMR90_upregulated","chr"=chrL,"First"=FL,"Second"=SL))
            }else if(gene[kk,]$log2.fold_change<=-0){# beacuse first sample is HUVEC,but testing from gene file log2(value_2/value_1) has taken
              gene1=rbind(gene1,cbind(gene[kk,],"regulation"="IMR90_downregulated","chr"=chrL,"First"=FL,"Second"=SL))
            }
          }else{
            gene1=rbind(gene1,cbind(gene[kk,],"regulation"="insignificant","chr"=chrL,"First"=FL,"Second"=SL))
          }
        }
      }
      saveRDS(gene1,paste(folder2,"/",dataName1,".rds",sep=""))
      
      ################################################################################################################
      
      if(!file.exists(paste(folder,"/","CM_",res,"_",pval,normMethod,"Fisher.rds",sep=""))){
        for(id in id_list){
          genesL=gene1[gene1$chr==id,]
          stripesL=stripes1[stripes1$chrom==id,]
          if(dim(genesL)[1]>0 & dim(stripesL)[1]>0){
            for(ii in 1:dim(genesL)[1]){
              GF=min(genesL$First[ii],genesL$Second[ii])
              GL=max(genesL$First[ii],genesL$Second[ii])
              for(jj in 1:dim(stripesL)[1]){
                print(paste(id,dim(genesL)[1],ii,dim(stripesL)[1],jj,sep="_"))
                #if(length(intersect(stripesL$firstLoc[[jj]]:stripesL$secondLoc[[jj]],genesL$First[ii] :genesL$Second[ii]))>0){
                if((stripesL$firstLoc[[jj]]>=GF & stripesL$firstLoc[[jj]]<=GL) | (stripesL$secondLoc[[jj]]>=GF & stripesL$secondLoc[[jj]]<=GL) ){
                  cm[genesL$regulation[ii],stripesL$sample1StripeType[[jj]]]=cm[genesL$regulation[ii],stripesL$sample1StripeType[[jj]]]+1
                }
              }
            }
          }
        }
        cm[cm==0]=1
        saveRDS(cm,paste(folder,"/","CM_",res,"_",pval,normMethod,".rds",sep=""))
        FET=fisher.test(cm,simulate.p.value=TRUE,B=1e7)$p.value
        saveRDS(FET,paste(folder,"/","CM_",res,"_",pval,normMethod,"Fisher.rds",sep=""))
      }
      cm=readRDS(paste(folder,"/","CM_",res,"_",pval,normMethod,".rds",sep=""))
      FET=readRDS(paste(folder,"/","CM_",res,"_",pval,normMethod,"Fisher.rds",sep=""))
      df=c()
      for(rr in rownames(cm)){
        for(cc in colnames(cm)){
          df=rbind(df,c(rr,cc,cm[rr,cc]))
        }
      }
      colnames(df)=c("data1","data2","values")
      df=data.frame(df)
      df$values=as.numeric(df$values)
      pdf(paste(folder,"/","CM_",res,"_",pval,normMethod,"confusionPlot.pdf",sep=""),height = 10, width = 10)
      print(ggplot(data =  df, mapping = aes(x = data1, y = data2)) +
              geom_tile(aes(fill = values), colour = "white") +
              geom_text(aes(label = values), vjust = 1) +
              scale_fill_gradient(low = "white", high = "steelblue")+
              scale_x_discrete(name="gene type", breaks=waiver(), labels=waiver(), limits=c("IMR90_upregulated","IMR90_downregulated","insignificant"))+
              scale_y_discrete(name="stripe type", breaks=waiver(), labels=waiver(), limits=c("insignificant","IMR90_weakened","IMR90_strengthened"))+
              theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size=15))+
              theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
      print(dev.off())
      pdf(paste(folder,"/","CM_",res,"_",pval,normMethod,"confusionPlot1.pdf",sep=""),height = 10, width = 10)
      print(ggplot(data =  df, mapping = aes(x = data1, y = data2)) +
              geom_tile(aes(fill = values), colour = "white") +
              geom_text(aes(label = values), vjust = 1) +
              scale_fill_gradient(low = "white", high = "steelblue")+
              ggtitle(paste("Fisher's exact test : ",FET,sep=""))+
              scale_x_discrete(name="gene type", breaks=waiver(), labels=waiver(), limits=c("IMR90_upregulated","IMR90_downregulated","insignificant"))+
              scale_y_discrete(name="stripe type", breaks=waiver(), labels=waiver(), limits=c("insignificant","IMR90_weakened","IMR90_strengthened"))+
              theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size=15))+
              theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
      print(dev.off())
      pdf(paste(folder,"/../","CM_",res,"_",pval,normMethod,"confusionPlot_6_1.pdf",sep=""),height = 10, width = 10)
      print(ggplot(data =  df, mapping = aes(x = data1, y = data2)) +
              geom_tile(aes(fill = values), colour = "white") +
              geom_text(aes(label = values), vjust = 1) +
              scale_fill_gradient(low = "white", high = "steelblue")+
              ggtitle(paste("Fisher's exact test : ",FET,sep=""))+
              scale_x_discrete(name="gene type", breaks=waiver(), labels=waiver(), limits=c("IMR90_upregulated","IMR90_downregulated","insignificant"))+
              scale_y_discrete(name="stripe type", breaks=waiver(), labels=waiver(), limits=c("insignificant","IMR90_weakened","IMR90_strengthened"))+
              theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size=15))+
              theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
      print(dev.off())
      ########################################################## another way
      ########################################## IMR90
      IMR90=matrix(0,ncol=2,nrow=2)
      rownames(IMR90)=c("IMR90_upregulated","IMR90_not_upregulated")
      colnames(IMR90)=c("IMR90_strengthened","IMR90_not_strengthened")
      IMR90=data.frame(IMR90)
      IMR90["IMR90_upregulated","IMR90_strengthened"]=cm["IMR90_upregulated","IMR90_strengthened"]
      IMR90["IMR90_upregulated","IMR90_not_strengthened"]=sum(cm["IMR90_upregulated",c("IMR90_weakened","insignificant")])
      IMR90["IMR90_not_upregulated","IMR90_strengthened"]=sum(cm[c("IMR90_downregulated","insignificant"),"IMR90_strengthened"])
      IMR90["IMR90_not_upregulated","IMR90_not_strengthened"]=(dim(gene)[1]+dim(stripes1)[1])-
      (cm["IMR90_upregulated","IMR90_strengthened"]+sum(cm["IMR90_upregulated",c("IMR90_weakened","insignificant")])+sum(cm[c("IMR90_downregulated","insignificant"),"IMR90_strengthened"]))
      FET=fisher.test(IMR90,simulate.p.value=TRUE,B=1e7,alternative = "greater")$p.value
      IMR90_FET=FET
      saveRDS(FET,paste(folder,"/","CM_",res,"_",pval,normMethod,"Fisher_IMR90.rds",sep=""))
      saveRDS(IMR90,paste(folder,"/","CM_",res,"_",pval,normMethod,"IMR90.rds",sep=""))
      df=c()
      for(rr in rownames(IMR90)){
        for(cc in colnames(IMR90)){
          df=rbind(df,c(rr,cc,IMR90[rr,cc]))
        }
      }
      colnames(df)=c("data1","data2","values")
      df=as.data.frame(df)
      df$values=as.numeric(df$values)
      pdf(paste(folder,"/","CM_",res,"_",pval,normMethod,"confusionPlot_IMR90.pdf",sep=""),height = 10, width = 10)
      print(ggplot(data =  df, mapping = aes(x = data1, y = data2)) +
              geom_tile(aes(fill = values), colour = "white") +
              geom_text(aes(label = values), vjust = 1) +
              scale_fill_gradient(low = "white", high = "steelblue")+
              ggtitle(paste("Fisher's exact test : ",FET,sep=""))+
              scale_y_discrete(name="gene type", breaks=waiver(), labels=waiver(), limits=c("IMR90_not_strengthened","IMR90_strengthened"))+
              scale_x_discrete(name="stripe type", breaks=waiver(), labels=waiver(), limits=c("IMR90_upregulated","IMR90_not_upregulated"))+
              theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size=15))+
              theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
      print(dev.off())
      pdf(paste(folder,"/../","CM_",res,"_",pval,normMethod,"confusionPlot_6_IMR90.pdf",sep=""),height = 10, width = 10)
      print(ggplot(data =  df, mapping = aes(x = data1, y = data2)) +
              geom_tile(aes(fill = values), colour = "white") +
              geom_text(aes(label = values), vjust = 1) +
              scale_fill_gradient(low = "white", high = "steelblue")+
              ggtitle(paste("Fisher's exact test : ",FET,sep=""))+
              scale_y_discrete(name="gene type", breaks=waiver(), labels=waiver(), limits=c("IMR90_not_strengthened","IMR90_strengthened"))+
              scale_x_discrete(name="stripe type", breaks=waiver(), labels=waiver(), limits=c("IMR90_upregulated","IMR90_not_upregulated"))+
              theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size=15))+
              theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
      print(dev.off())
      
      
      ########################################## HUVEC
      HUVEC=matrix(0,ncol=2,nrow=2)
      rownames(HUVEC)=c("HUVEC_upregulated","HUVEC_not_upregulated")
      colnames(HUVEC)=c("HUVEC_strengthened","HUVEC_not_strengthened")
      HUVEC=as.data.frame(HUVEC)
      HUVEC["HUVEC_upregulated","HUVEC_strengthened"]=cm["IMR90_downregulated","IMR90_weakened"]
      HUVEC["HUVEC_upregulated","HUVEC_not_strengthened"]=sum(cm["IMR90_downregulated",c("IMR90_strengthened","insignificant")])
      HUVEC["HUVEC_not_upregulated","HUVEC_strengthened"]=sum(cm[c("IMR90_upregulated","insignificant"),"IMR90_weakened"])
      HUVEC["HUVEC_not_upregulated","HUVEC_not_strengthened"]=(dim(gene)[1]+dim(stripes1)[1])-
      (cm["IMR90_downregulated","IMR90_weakened"]+sum(cm["IMR90_downregulated",c("IMR90_strengthened","insignificant")])+sum(cm[c("IMR90_upregulated","insignificant"),"IMR90_weakened"]))
      FET=fisher.test(HUVEC,simulate.p.value=TRUE,B=1e7,alternative = "greater")$p.value
      HUVEC_FET=FET
      saveRDS(FET,paste(folder,"/","CM_",res,"_",pval,normMethod,"Fisher_HUVEC.rds",sep=""))
      saveRDS(HUVEC,paste(folder,"/","CM_",res,"_",pval,normMethod,"HUVEC.rds",sep=""))
      df=c()
      for(rr in rownames(HUVEC)){
        for(cc in colnames(HUVEC)){
          df=rbind(df,c(rr,cc,HUVEC[rr,cc]))
        }
      }
      colnames(df)=c("data1","data2","values")
      df=as.data.frame(df)
      df$values=as.numeric(df$values)
      pdf(paste(folder,"/","CM_",res,"_",pval,normMethod,"confusionPlot_HUVEC.pdf",sep=""),height = 10, width = 10)
      print(ggplot(data =  df, mapping = aes(x = data1, y = data2)) +
              geom_tile(aes(fill = values), colour = "white") +
              geom_text(aes(label = values), vjust = 1) +
              scale_fill_gradient(low = "white", high = "steelblue")+
              ggtitle(paste("Fisher's exact test : ",FET,sep=""))+
              scale_y_discrete(name="gene type", breaks=waiver(), labels=waiver(), limits=c("HUVEC_not_strengthened","HUVEC_strengthened"))+
              scale_x_discrete(name="stripe type", breaks=waiver(), labels=waiver(), limits=c("HUVEC_upregulated","HUVEC_not_upregulated"))+
              theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size=15))+
              theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
      print(dev.off())
      pdf(paste(folder,"/../","CM_",res,"_",pval,normMethod,"confusionPlot_6_HUVEC.pdf",sep=""),height = 10, width = 10)
      print(ggplot(data =  df, mapping = aes(x = data1, y = data2)) +
              geom_tile(aes(fill = values), colour = "white") +
              geom_text(aes(label = values), vjust = 1) +
              scale_fill_gradient(low = "white", high = "steelblue")+
              ggtitle(paste("Fisher's exact test : ",FET,sep=""))+
              scale_y_discrete(name="gene type", breaks=waiver(), labels=waiver(), limits=c("HUVEC_not_strengthened","HUVEC_strengthened"))+
              scale_x_discrete(name="stripe type", breaks=waiver(), labels=waiver(), limits=c("HUVEC_upregulated","HUVEC_not_upregulated"))+
              theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size=15))+
              theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
      print(dev.off())
      ########################################## combined bar plot
      # # random control
      # bar_df=data.frame("percentage"=c(HUVEC["HUVEC_upregulated","HUVEC_strengthened"]/(HUVEC["HUVEC_upregulated","HUVEC_strengthened"]+HUVEC["HUVEC_not_upregulated","HUVEC_strengthened"]),
      #                                  (HUVEC["HUVEC_upregulated","HUVEC_strengthened"]+HUVEC["HUVEC_upregulated","HUVEC_not_strengthened"])/sum(HUVEC),
      #                                  IMR90["IMR90_upregulated","IMR90_strengthened"]/(IMR90["IMR90_upregulated","IMR90_strengthened"]+IMR90["IMR90_not_upregulated","IMR90_strengthened"]),
      #                                  (IMR90["IMR90_upregulated","IMR90_strengthened"]+IMR90["IMR90_upregulated","IMR90_not_strengthened"])/sum(IMR90))*100,
      #                   "stripeType"=c("HUVEC strengthened stripe","Random control of HUVEC","IMR90 strengthened stripe","Random control of IMR90"),
      #                   "geneType"=c("Up-regulated in HUVEC","Up-regulated in HUVEC","Up-regulated in IMR90","Up-regulated in IMR90"))
      
      # expected percentage
      bar_df=data.frame("percentage"=c(HUVEC["HUVEC_upregulated","HUVEC_strengthened"]/(HUVEC["HUVEC_upregulated","HUVEC_strengthened"]+HUVEC["HUVEC_not_upregulated","HUVEC_strengthened"]),
                                       ((HUVEC["HUVEC_upregulated","HUVEC_strengthened"]+HUVEC["HUVEC_not_upregulated","HUVEC_strengthened"])*
                                                     (HUVEC["HUVEC_upregulated","HUVEC_strengthened"]+HUVEC["HUVEC_upregulated","HUVEC_not_strengthened"])/sum(HUVEC))/
                                         (HUVEC["HUVEC_upregulated","HUVEC_strengthened"]+HUVEC["HUVEC_not_upregulated","HUVEC_strengthened"]),
                                       IMR90["IMR90_upregulated","IMR90_strengthened"]/(IMR90["IMR90_upregulated","IMR90_strengthened"]+IMR90["IMR90_not_upregulated","IMR90_strengthened"]),
                                       ((IMR90["IMR90_upregulated","IMR90_strengthened"]+IMR90["IMR90_not_upregulated","IMR90_strengthened"])*
                                                     (IMR90["IMR90_upregulated","IMR90_strengthened"]+IMR90["IMR90_upregulated","IMR90_not_strengthened"])/sum(IMR90))/
                                         (IMR90["IMR90_upregulated","IMR90_strengthened"]+IMR90["IMR90_not_upregulated","IMR90_strengthened"])
      )*100,
      "stripeType"=c("HUVEC strengthened stripe","Random control of HUVEC","IMR90 strengthened stripe","Random control of IMR90"),
      "geneType"=c("Up-regulated in HUVEC","Up-regulated in HUVEC","Up-regulated in IMR90","Up-regulated in IMR90"))

      bar_pval_df=rbind(bar_pval_df,c(normMethod,res,pval,"HUVEC",bar_df[1,"percentage"],bar_df[2,"percentage"],round(bar_df[1,"percentage"]/bar_df[2,"percentage"],2)))
      bar_pval_df=rbind(bar_pval_df,c(normMethod,res,pval,"IMR90",bar_df[3,"percentage"],bar_df[4,"percentage"],round(bar_df[3,"percentage"]/bar_df[4,"percentage"],2)))
      saveRDS(bar_pval_df,paste(folder,"/../confusionPlot_6_varingNorm_HUVEC_IMR90.rds",sep=""))
      
      temp1=strsplit(as.character(HUVEC_FET),"e")
      temp2=strsplit(temp1[[1]][1],".")
      HUVEC_FET=paste(substr(temp1[[1]][1],1,3),"e",temp1[[1]][2],sep="")
      temp1=strsplit(as.character(IMR90_FET),"e")
      temp2=strsplit(temp1[[1]][1],".")
      IMR90_FET=paste(substr(temp1[[1]][1],1,3),"e",temp1[[1]][2],sep="")
      
      
      # pdf(paste(folder,"/../","CM_",res,"_",pval,normMethod,"confusionPlot_6_HUVEC_IMR90.pdf",sep=""),height = 10, width = 10)
      print(ggplot(bar_df, aes(x=stripeType, y=percentage, fill=geneType)) + geom_bar(stat="identity", position=position_dodge())+
              theme_classic()+ 
              # ggtitle("CTCF binding motifs orientation at stripe")+ #xlab("Stripe's Types") + ylab("Orientation Percentage")+ 
              scale_x_discrete(name="Type", breaks=waiver(), labels=waiver(), limits=bar_df$stripeType)+
              scale_y_continuous(name="Percentage", breaks=waiver(), labels=waiver(), limits=NULL)+
              theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size=10,family = "Helvetica"))+
              annotate(geom="text",x=c(1.5, 3.5),y=max(bar_df$percentage)+1,label=c(paste("p = ",HUVEC_FET,sep=""),paste("p = ",IMR90_FET,sep="")),color="black",size=2)+
              theme(legend.position="top",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)))
      ggsave(filename=paste(folder,"/../","CM_",res,"_",pval,normMethod,"confusionPlot_6_HUVEC_IMR90.svg",sep=""),dpi = 300,width=2,height=4)
      print(dev.off())
      
    }  
  }  
}

library(plyr)
colnames(bar_pval_df)<-bar_pval_df_c
bar_pval_df$fold=as.numeric(bar_pval_df$fold)
bar_pval_df$data=factor(bar_pval_df$data,c("HUVEC","IMR90"))
bar_pval_df1=bar_pval_df
# bar_pval_df1 <- arrange(bar_pval_df1,data,normMethod)
bar_pval_df1 <- ddply(bar_pval_df1, "normMethod",transform, label_ypos=cumsum(fold))
head(bar_pval_df1)
print(ggplot(bar_pval_df1, aes(x=normMethod, y=fold, fill=data)) + geom_bar(stat="identity")+
        theme_classic()+ 
        # ggtitle("CTCF binding motifs orientation at stripe")+ #xlab("Stripe's Types") + ylab("Orientation Percentage")+ 
        scale_x_discrete(name="Normalization methods", breaks=waiver(), labels=waiver(), limits=NULL)+
        scale_y_continuous(name="Fold change", breaks=waiver(), labels=waiver(), limits=NULL)+
        theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size=10,family = "Helvetica"))+
        # geom_text(aes(y=label_ypos, label=fold), vjust=1.6, color="white", size=3.5)+
        theme(legend.position="top",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+facet_grid(~data))
ggsave(filename=paste(folder,"/../","barplot_","confusionPlot_6_varingNorm_HUVEC_IMR90.svg",sep=""),dpi = 300,width=3,height=4)
print(dev.off())
###########################################

###########################################

###########################################

###########################################
