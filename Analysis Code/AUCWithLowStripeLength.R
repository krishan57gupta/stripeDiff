library(sigmoid)
library(caret)
library(pROC)
library(ggplot2)
library(export)
library(Hmisc)
library(ggpubr)
library("cowplot")
library("gridExtra")
N<-function(x,minv=0,maxv=1){
  (((x-min(x))/(max(x)-min(x)))*(maxv-minv))+minv
}
folder="~/Revision/2_2_181"

############################################################
windows=c(seq(20,90,10))
score4=c()
rocScoreAll=list()
rocPlotAll=list()
for(w in windows){
  rocScore2=list()
  data=read.csv(paste(folder,"/result_10000_10-110L_",w,".txt",sep=""), sep = "\t")
  val=round(auc(data$label, data$score),2)
  score4=c(score4,val)
  rocobj1 <- smooth(roc(data$label, data$score))
  rocScoreAll[[paste("Window Size: 6MB, Stripe Length : ",w/100,sep="")]]=rocobj1
  rocScore2[[paste("AUC : ",val,sep="")]]=rocobj1
 
  data=read.csv(paste(folder,"/result_10000_10-190L_",w,".txt",sep=""), sep = "\t")
  val=round(auc(data$label, data$score),3)
  score4=c(score4,val)
  rocobj2 <- smooth(roc(data$label, data$score))
  rocScoreAll[[paste("Window Size: 20MB, Stripe Length : ",w/100,sep="")]]=rocobj2
  rocScore2[[paste("AUC : ",val,sep="")]]=rocobj2
 
  rocPlotAll[[paste("Stripe Length : ",w/100," MB",sep="")]]<- ggroc(rocScore2)+
    theme(text = element_text(size=10,family = "Helvetica"))+
    ggtitle(paste(" Stripe Length : ",w/100," MB",sep=""))+theme_classic()+
    theme(legend.position = "none")
}

plotSW <- do.call(grid.arrange, list(grobs=rocPlotAll, ncol=3, width=3, label.x="specificity",label.y="sensitivity", common.legend=TRUE))
plot(plotSW)
df=data.frame("AUC"=score4,"length"=as.factor(sort(rep(windows/100,2))),
              label=rep(c("Window size = 6MB","Window size = 20MB"),length(windows)))#[sort(c(c(2,4,6,8)*2,c(2,4,6,8)*2-1)),]
# pdf(paste(folder,"/AUCwithStripeLength_10KB.pdf",sep=""),width=10,height = 10)
p=ggplot(df, aes(x=length, y=AUC,color=label)) + geom_point(size = 1)+theme_classic()+ xlab("Stripe size in (MB)") + ylab("AUC")+
  scale_x_discrete(name="Stripe size in (MB)", breaks=waiver(), labels=waiver(), limits=as.factor(windows/100))+
  scale_y_continuous(name="AUC", breaks=waiver(), labels=waiver(), limits=c(0.6,1))+
  theme(text = element_text(size=10,family = "Helvetica"))+
  theme(legend.position="top")
update_labels(p, list(colour=""))
ggsave(filename=paste(folder,"/AUCwithStripeLength_10KB.pdf",sep=""),dpi = 300,width=5,height=5)
dev.off()


plot(plotSW)
ggsave(filename=paste(folder,"/AUROCwithStripeLength_10KB.pdf",sep=""),dpi = 300,width=5,height=5)
dev.off()


