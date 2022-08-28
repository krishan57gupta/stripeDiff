folder="~/Revision/1_12/"
library(export)
library(Hmisc)
library(ggplot2)

############### correlation
folder="~/Revision/1_12/"
if(file.exists(paste(folder,"/IMR90_HUVEC_10kb_121_123.txt",sep=""))){
  df=read.table(paste(folder,"/IMR90_HUVEC_10kb_121_123.txt",sep=""))
  start=dim(df)[1]+1
}else{
  df=c()
  start=length(df)+1
}
data_IMR90_HUVEC_chr_10kb <- readRDS(paste("~/Revision/2_4_1/data_VC1/data_IMR90_HUVEC_chr8_10kb.rds",sep=""))[[1]]
dim(data_IMR90_HUVEC_chr_10kb)
data2=data_IMR90_HUVEC_chr_10kb[12101:12300,12101:12300]
dim(data2)
# library(pheatmap)
# pheatmap::pheatmap(data,cluster_rows = F,cluster_cols = F)

source('~/Desktop/stripeDiff-master/stripeDiff_Krishan/Revision/1_12/simulation_stripe.R')
score=c()
upLimit=100
for(kk in start:10000){
  ntad=5
  lt=100
  ft=10
  bin=10000
  outputPath=paste('~/Revision/1_12/',sep="")
  generateSimulationData(4,outputPath)
  data1=read.table("~/Revision/1_12/simulation4_pair1.txt")
  data1[data1>upLimit]=upLimit
  data2[data2>upLimit]=upLimit
  dim(data1)
  k=min(dim(data2)[1],dim(data1)[1])
  data2=data2[1:k,1:k]
  data1=as.numeric(unlist(data1[1:k,1:k]))
  score=cor(c(data2),as.numeric(unlist(data1)))
  print(paste(kk,score,mean(data2),mean(data1),sep="_"))
  df=rbind(df,c(score,"chr8: 121-123 mb"))
  write.table(df,paste(folder,"/IMR90_HUVEC_10kb_121_123.txt",sep=""))
}
df=read.table(paste(folder,"IMR90_HUVEC_10kb_121_123.txt",sep=""))[1:100,]
colnames(df)=c("pearson","IMR90")

plot=ggplot(df, aes(x=IMR90, y=pearson, fill=IMR90)) +
  geom_boxplot()+
  theme_classic()+ ggtitle(paste("",sep=""))+
  scale_x_discrete(name="IMR90", breaks=waiver(), labels=waiver(), limits=NULL)+
  scale_y_continuous(name="Pearson correlation", breaks=waiver(), labels=waiver(), limits=c(0.7,.85))+
  theme(plot.title = element_text(hjust = 0.5))+ 
  theme(text = element_text(size=10,family = "Helvetica"))+
  theme(legend.position="top",axis.text.x = element_text(angle = 0, vjust = .5, hjust=.5))# +coord_flip()
plot
ggsave(filename=paste(folder,"pearson_IMR90_HUVEC_10kb_121_123",".pdf",sep=""),dpi = 300,width=5,height=5)
dev.off()



