library(sigmoid)
library(caret)
library(pROC)
library(export)
library(Hmisc)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(cowplot)


N<-function(x,minv=0,maxv=1){
  (((x-min(x))/(max(x)-min(x)))*(maxv-minv))+minv
}
folder="~/Revision/2_1_2"
############################################################################################


################################# ROC plot
####################### comparing with edgeR
bin=10000
ft=10
lt=110
# windows=c(seq(0,50,5),seq(2,8,2))
windows=factor(c(0,2,4,5,6,8,10,15,20,25,30,35,40,45,50))
# windows=c(seq(0,50,5))
for(noise in c("GD","RL")[1:2]){
  score4=c()
  rocScoreAll=list()
  rocPlotAll=list()
  for(i in windows){
    rocScore2=list()
    w=300
    data=read.csv(paste(folder,'/result_',bin,'_',ft,'-',lt,'_',noise,'_',i,'_',w,'_DF_DS_new7.txt',sep=""), sep = "\t")
    val=round(auc(data$label, data$score),2)
    score4=c(score4,val)
    rocobj1 <- smooth(roc(data$label, data$score))
    rocScoreAll[[paste("Method : stripeDiff, ",noise," Noise : ",i,sep="")]]=rocobj1
    rocScore2[[paste("AUC : ",val,sep="")]]=rocobj1
    
    data=read.csv(paste(folder,'/result_',bin,'_',ft,'-',lt,'_',noise,'_',i,'_',w,'_DF_ER_new7.txt',sep=""), sep = "\t")
    val=round(auc(data$label, data$score),2)
    score4=c(score4,val)
    rocobj2 <- smooth(roc(data$label, data$score))
    rocScoreAll[[paste("Method : edgeR, ",noise," Noise : ",i,sep="")]]=rocobj2
    rocScore2[[paste("AUC : ",val,sep="")]]=rocobj2
    rocPlotAll[[paste(noise," Noise : ",i,"%",sep="")]]<- ggroc(rocScore2)+
      theme(text = element_text(size=10,family = "Helvetica"))+
      ggtitle(paste(noise," Noise : ",i,"%",sep=""))+theme_classic()+
      theme(legend.position = "none")
    # rocPlotAll[[paste(noise," Noise : ",i,"%",sep="")]]<- update_labels(rocPlotAll[[paste(noise," Noise : ",i,"%",sep="")]], list(colour=""))
  }

  plotSW <- do.call(grid.arrange, list(grobs=rocPlotAll, ncol=3, width=6, label.x="specificity",label.y="sensitivity", common.legend = TRUE))
  ggsave(filename=paste(folder,"/AUROCwithnoisesEdgeR_DS2_",noise,".pdf",sep=""),dpi = 300,width=5,height=5)
    print(dev.off())
  # plot(plotSW)
  
  df=data.frame("AUC"=score4,"length"=as.factor(sort(rep(windows,2))),
                "methods"=rep(c("DiffStripe","EdgeR"),length(windows)))#[sort(c(c(2,4,6,8)*2,c(2,4,6,8)*2-1)),]
  df$methods=factor(df$methods)
  # pdf(paste(folder,"/AUCwithnoisesEdgeR_DS2_",noise,".pdf",sep=""),width = 10,height =10 )
  print(ggplot(data=df, aes(x=length, y=AUC, group=methods)) +geom_line(aes(color=methods),size=.2)+geom_point(size=1,aes(color=methods))+theme_classic()+ggtitle(paste("Window size : 6MB and noise : ",noise,sep=""))+
          theme(plot.title = element_text(hjust = 0.5))+
          scale_x_discrete(name="Noise in Percentage", breaks=waiver(), labels=waiver(), limits=windows)+
          scale_y_continuous(name="AUC", breaks=waiver(), labels=waiver(), limits=c(0.5,1))+
          theme(text = element_text(size=10,family = "Helvetica"))+
          theme(legend.position="top"))
  ggsave(filename=paste(folder,"/AUCwithnoisesEdgeR_DS2_",noise,".pdf",sep=""),dpi = 300,width=5,height=5)
   print(dev.off())
}


###################################################### just comparing distances from actual
bin=10000
ft=10
lt=110
windows=factor(c(0))

for(noise in c("GD")){
  score4=c()
  scoreAll=c()
  for(i in windows){
    w=300
    data=read.csv(paste(folder,'/result_',bin,'_',ft,'-',lt,'_',noise,'_',i,'_',w,'_DF_DSF_new10.txt',sep=""), sep = "\t")
    data$score=log2(data$score/bin+1)
    temp=cbind(data,"method"=rep('diffStripe',dim(data)[1]),"noisePercentage"=factor(rep(i,dim(data)[1])))
    score4=rbind(score4,temp)
    
    # temp=cbind("label"=1,"median1"=median(data[data$label==1,2]),"mad1"=mad(data[data$label==1,2]),"method"='diffStripe',"noisePercentage"=i)
    # scoreAll=rbind(scoreAll,temp)
    temp=cbind("label"=0,"median1"=median(data[data$label==0,2]),"mad1"=mad(data[data$label==0,2]),"method"='diffStripe',"noisePercentage"=i)
    scoreAll=rbind(scoreAll,temp)
    
    data=read.csv(paste(folder,'/result_',bin,'_',ft,'-',lt,'_',noise,'_',i,'_',w,'_DF_ZM_new10.txt',sep=""), sep = "\t")
    data$score=log2(data$score/bin+1)
    temp=cbind(data,"method"=rep('Zebra',dim(data)[1]),"noisePercentage"=factor(rep(i,dim(data)[1])))
    score4=rbind(score4,temp)
    
    # temp=cbind("label"=1,"median1"=median(data[data$label==1,2]),"mad1"=mad(data[data$label==1,2]),"method"='Zebra',"noisePercentage"=i)
    # scoreAll=rbind(scoreAll,temp)
    temp=cbind("label"=0,"median1"=median(data[data$label==0,2]),"mad1"=mad(data[data$label==0,2]),"method"='Zebra',"noisePercentage"=i)
    scoreAll=rbind(scoreAll,temp)
  }
 score4$noisePercentage=factor(score4$noisePercentage,levels = sort(windows))
 scoreAll=data.frame(scoreAll)
 scoreAll$noisePercentage=factor(scoreAll$noisePercentage,levels = sort(windows))
 scoreAll$median1=as.numeric(scoreAll$median1)
 scoreAll$mad1=as.numeric(scoreAll$mad1)


    score4_1=score4[score4$noisePercentage==0,]
    scoreAll_1=scoreAll[scoreAll$noisePercentage==0,]
    # pdf(paste(folder,"/DistanceWithnoisesZebra_DSF2_sim.pdf",sep=""),width = 10,height =10 )
    print(ggplot(data=score4_1, aes(x=method, y=score, fill=method)) +geom_boxplot(outlier.colour="white", outlier.shape=16,outlier.size=.01)+
            theme_classic()+ggtitle(paste("Window size : 6MB",sep=""))+
            theme(plot.title = element_text(hjust = 0.5))+
            scale_x_discrete(name="method", breaks=waiver(), labels=waiver(), limits=NULL)+
            scale_y_continuous(name="Distance in 10kb", breaks=waiver(), labels=waiver(), limits=c(0,.2))+
            theme(text = element_text(size=10,family = "Helvetica"))+
            theme(legend.position="top",axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=.5))+
            stat_compare_means(label.x = 1.5))
            ggsave(filename=paste(folder,"/DistanceWithnoisesZebra_DSF2_sim.pdf",sep=""),dpi = 300,width=5,height=5)
    print(dev.off())
  
}

############################################################################################








