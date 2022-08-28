library(export)
library(Hmisc)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(gridGraphics)
stripeDiffDF=read.table(paste("~/Revision/1_13/stripeDiffDF.txt",sep=""))

# dff1=stripeDiffDF[stripeDiffDF$resolution=="10kb",]
dff1=stripeDiffDF

dff1$resolution=factor(dff1$resolution,levels = c("100kb","50kb","25kb","10kb"))

limit=5
pointSize=.1
x=scale_x_continuous(name="x", breaks=waiver(), labels=waiver(), limits=c(-limit,limit))
y=scale_y_continuous(name="y", breaks=waiver(), labels=waiver(), limits=c(-limit,limit))

A=ggplot(dff1, aes(x=down.tscore)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white",size=pointSize)+
  geom_density(alpha=.2, fill="red",size=pointSize) +
  # geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1))+
  # scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  theme_classic()+ ggtitle(paste("",sep=""))+ # xlab("hg19's chromosomes") + ylab("density")+ 
  x+
  theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size=10,family = "Helvetica"))+
  theme(legend.position="top",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))# +coord_flip()
A=A+facet_grid(~resolution)
A

B=ggplot(dff1, aes(sample=down.tscore)) + 
  stat_qq(size=pointSize) + stat_qq_line(size=pointSize)+
  # geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1))+
  # scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  theme_classic()+ ggtitle(paste("",sep=""))+#  xlab("hg19's chromosomes") + ylab("density")+ 
  y+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size=10,family = "Helvetica"))+
  theme(legend.position="top",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))# +coord_flip()
B=B+facet_grid(~resolution)
B

##########################################################################################################

dff3=stripeDiffDF
for(kk in c("100kb","50kb","25kb","10kb")){
  index=dff3$resolution==kk
  dff2=dff3[index,]
  u1 = dff2$upPeak.sample1
  u2 = dff2$upPeak.sample2
  m = min(u1, u2)
  if(m<0){
    u1 = u1-1.1*m
    u2 = u2-1.1*m
  }
  x1 = u1/u2
  m1 = mean(x1)
  sd1 = sd(x1)
  
  d1 = dff2$downPeak.sample1
  d2 = dff2$downPeak.sample2
  m = max(d1,d2)
  if(m>0){
    d1 = d1-1.1*m
    d2 = d2-1.1*m
  }
  x2 = -d1/(-d2)
  m2 = mean(x2)
  sd2 = sd(x2)
  
  # x1.norm = (x1-m1)/sd1
  # x2.norm = (x2-m2)/sd2
  # chiq = x1.norm^2+x2.norm^2
  # stripe$chiq = chiq
  # stripe$p.chip = 1-pchisq(chiq,1)
  
  dff3$up.tscore[index] = (x1-m1)/(sd1)
  # stripe$p.up = 2*pt(-abs(stripe$up.tscore),df=nrow(stripe)-1)
  dff3$down.tscore[index]  = (x2-m2)/(sd2)
  # stripe$p.down = 2*pt(-abs(stripe$down.tscore),df=nrow(stripe)-1)
  ks_up=round(ks.test(dff3$up.tscore[index],rnorm(100,mean = 0,sd=1))$p.value,2)
  ks_down=round(ks.test(dff3$up.tscore[index],rnorm(100,mean = 0,sd=1))$p.value,2)
  print(paste(kk,ks_up,ks_down,sep="_"))
}

dff3$resolution=factor(dff3$resolution,levels = c("100kb","50kb","25kb","10kb"))


C=ggplot(dff3, aes(x=down.tscore)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white",size=pointSize)+
  geom_density(alpha=.2, fill="red",size=pointSize) +
  # geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1))+
  # scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  theme_classic()+ ggtitle(paste("",sep=""))+ # xlab("hg19's chromosomes") + ylab("density")+ 
  x+
  theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size=10,family = "Helvetica"))+
  theme(legend.position="top",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))# +coord_flip()
C=C+facet_grid(~resolution)
C

D=ggplot(dff3, aes(sample=down.tscore)) + 
  stat_qq(size=pointSize) + stat_qq_line(size=pointSize)+
  # geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1))+
  # scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  theme_classic()+ ggtitle(paste("",sep=""))+#  xlab("hg19's chromosomes") + ylab("density")+ 
  y+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size=10,family = "Helvetica"))+
  theme(legend.position="top",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))# +coord_flip()
D=D+facet_grid(~resolution)
D

ABCD <- ggarrange(A,B,C,D, labels = c("A","B","C","D"),ncol = 2, nrow = 2)
ABCD
ggsave(filename=paste("~/Desktop/stripeDiff-master/stripeDiff_Krishan/Revision/1_11/normal",limit,".pdf",sep=""),dpi = 300,width=5,height=5)
dev.off()