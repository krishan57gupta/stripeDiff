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
folder="~/Revision/2_2_10_Latest"
dataName="IMR90"
seqD=factor(seq(50,250,50))


####################### comparing
windows=c(seq(20,300,20))
score4=c()
rocScoreAll=list()
rocPlotAll=list()
for(w in windows){
  rocScore4=list()
  for(seq in seqD){
    data=read.csv(paste(folder,"/sim",dataName,"_",seq,"_",w,'.txt',sep=""), sep = "\t")
    val=round(auc(data$label, data$score),2)
    score4=c(score4,val)
    rocobj <-smooth(
      roc(data$label, data$score)
    )
    rocScoreAll[[paste("Depth Size: ",seqD[seq],"MB", " Stripe Length : ",w/100,sep="")]]=rocobj
    rocScore4[[paste("AUC : ",val,sep="")]]=rocobj
    
  }
  rocPlotAll[[paste("Stripe Length : ",w/100,sep="")]]<- ggroc(rocScore4,size = .2)+
    theme(text = element_text(size=10,family = "Helvetica"))+
    ggtitle(paste(" Stripe Length : ",w/100," MB",sep=""))+theme_classic()+
    theme(legend.position = "none")
}
plotSW <- do.call(grid.arrange, list(grobs=rocPlotAll, ncol=3, width=6, label.x="specificity",label.y="sensitivity", common.legend = TRUE))
df=data.frame("AUC"=score4,"length"=as.factor(sort(rep(windows/100,length(seqD)))),
              depth=rep(c(seqD),length(windows)))#[sort(c(c(2,4,6,8)*2,c(2,4,6,8)*2-1)),]
p=ggplot(df, aes(x=length, y=AUC, group=depth)) + 
  geom_point(aes(color=depth),size = 1)+
  geom_line(aes(color=depth),size=.2)+  # ggtitle(paste("Depths in MB",sep=""))+
  theme_classic()+
  scale_x_discrete(name="Stripe size in MB", breaks=waiver(), labels=waiver(), limits=NULL)+
  scale_y_continuous(name="AUC", breaks=waiver(), labels=waiver(), limits=NULL)+
  theme(text = element_text(size=10,family = "Helvetica"))+
  theme(legend.position="top",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=.5))
  update_labels(p, list(colour="Depth in MB"))
  ggsave(filename=paste(folder,"/AUCwithStripeLength_Depth_10KB.pdf",sep=""),dpi = 300,width=5,height=5)
  graph2ppt(file=paste("~/Desktop/stripeDiff-master/stripeDiff_Krishan/Revision/diffStripe.pptx",sep=""), width=3, height=3, append=TRUE) 
dev.off()
# pdf(paste(folder,"/AUROCwithStripeLength_Depth_10KB.pdf",sep=""),width=10,height = 10)
plot(plotSW)
ggsave(filename=paste(folder,"/AUROCwithStripeLength_Depth_10KB.pdf",sep=""),dpi = 300,width=5,height=5)
dev.off()

