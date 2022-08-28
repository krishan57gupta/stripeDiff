library(sigmoid)
library(caret)
library(pROC)
library(export)
library(Hmisc)
library(ggpubr)
library("cowplot")
library("gridExtra")
N<-function(x,minv=0,maxv=1){
  (((x-min(x))/(max(x)-min(x)))*(maxv-minv))+minv
}
folder="~/Revision/2_1_1"


windows=c(seq(20,90,10))
score4=c()
rocScoreAll=list()
rocPlotAll=list()
for(w in windows){
  rocScore2=list()
  data=read.csv(paste(folder,"/result_10000_10-110L_",w,"_DF.txt",sep=""), sep = "\t")
  val=round(auc(data$label, data$score),2)
  score4=c(score4,val)
  rocobj1 <- smooth(roc(data$label, data$score))
  rocScoreAll[[paste("Method : stripeDiff, "," stripe size : ",w/100,sep="")]]=rocobj1
  rocScore2[[paste("AUC : ",val,sep="")]]=rocobj1
  
  data=read.csv(paste(folder,"/result_10000_10-110L_",w,"_FC.txt",sep=""), sep = "\t")
  val=round(auc(data$label, data$score),2)
  score4=c(score4,val)
  rocobj2 <- smooth(roc(data$label, data$score))
  rocScoreAll[[paste("Method : logFoldChange, "," stripe size : ",w/100,sep="")]]=rocobj2
  rocScore2[[paste("AUC : ",val,sep="")]]=rocobj2
  rocPlotAll[[paste(" stripe size : ",w/100," MB",sep="")]]<- ggroc(rocScore2)+
    theme(text = element_text(size=10,family = "Helvetica"))+
    ggtitle(paste(" striep size : ",w/100," MB",sep=""))+theme_classic()+
    theme(legend.position = "none")
  # rocPlotAll[[paste(" stripe size : ",w/100," MB",sep="")]]<- update_labels(rocPlotAll[[paste(" stripe size : ",w/100," MB",sep="")]], list(colour=""))
}

plotSW <- do.call(grid.arrange, list(grobs=rocPlotAll, ncol=3, width=6, label.x="specificity",label.y="sensitivity", common.legend = TRUE))
# plot(plotSW)
ggsave(filename=paste(folder,"/AUROCwithStripeSizeLogFoldChabge_4.pdf",sep=""),dpi = 300,width=5,height=5)
graph2ppt(file=paste("~/Revision/diffStripe.pptx",sep=""), width=6, height=6, append=TRUE) 
dev.off()
df=data.frame("AUC"=score4,"length"=as.factor(sort(rep(windows/100,2))),
              methods=rep(c("DiffStripe","logFoldChange"),length(windows)))#[sort(c(c(2,4,6,8)*2,c(2,4,6,8)*2-1)),]
df$methods=factor(df$methods)
print(ggplot(data=df, aes(x=length, y=AUC, group=methods)) +geom_line(aes(color=methods))+geom_point(size=1,aes(color=methods))+
        theme_classic()+ggtitle(paste("Window size : 6MB",sep=""))+
        theme(plot.title = element_text(hjust = 0.5))+
        scale_x_discrete(name="Stripe size in (MB)", breaks=waiver(), labels=waiver(), limits=as.factor(windows/100))+
        scale_y_continuous(name="AUC", breaks=waiver(), labels=waiver(), limits=c(0.6,1))+
        theme(text = element_text(size=10,family = "Helvetica"))+
        theme(legend.position="top",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
ggsave(filename=paste(folder,"/AUCwithStripeSizeLogFoldChabge_4.pdf",sep=""),dpi = 300,width=5,height=5)
dev.off()




