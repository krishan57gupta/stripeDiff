library(readxl)
library(ggplot2)
library(ggpubr)
library(export)
library(Hmisc)
library(scales)
# scientific_10 <- function(x) {
#   parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
# }
scientific_10 <- function(x) {
  parse(text=gsub("e","%*%10^", scales::scientific_format()(x)))
}

folder="~/Revision/Figure_4E/"
HUVEC<-read_excel(paste(folder,"data strengthened in IMR90.xlsx",sep=""))
IMR90<-read_excel(paste(folder,"data strengthened in HUVEC.xlsx",sep=""))

HUVEC_strengthened=c()
for(i in 1:dim(HUVEC)[1]){
  for(j in 1:length(HUVEC[i,])){
    HUVEC_strengthened=rbind(HUVEC_strengthened,c(HUVEC[i,][j][[1]],c(strsplit(names(HUVEC[i,][j]),"_")[[1]]),i))
  }
}
colnames(HUVEC_strengthened)<-c("signal","protien","data","id")
HUVEC_strengthened=as.data.frame(HUVEC_strengthened)
HUVEC_strengthened$signal=as.numeric(HUVEC_strengthened$signal)

HUVEC_strengthened_Pol2=HUVEC_strengthened[HUVEC_strengthened$protien=="Pol2",]
plot<-ggplot(HUVEC_strengthened_Pol2, aes(x=data, y=signal, fill=data)) + geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=2)+
  theme_classic()+ ggtitle("Stripes strengthened in IMR90 relative to HUVEC cells")+
  scale_x_discrete(name="", breaks=waiver(), labels=waiver(), limits=NULL)+
  scale_y_continuous(name="Chip-seq signal", breaks=waiver(), labels=scientific_10, limits=NULL)+
  # scale_y_log10("Chip-seq signal",breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size=10,family = "Helvetica"))+
  stat_compare_means(method = "wilcox.test",comparisons = list(c("HUVEC","IMR90")),method.args=list(alternative = "greater"),paired=TRUE,
                     aes(label =  paste0(..method.., 'p =', ..p.format..)), label.y = 220)+
  theme(legend.position="top",axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=.5))
print(plot)
ggsave(filename=paste(folder,"HUVEC_Pol2.svg",sep=""),dpi = 300,width=5,height=5)
print(dev.off())

HUVEC_strengthened_H3K27ac=HUVEC_strengthened[HUVEC_strengthened$protien=="H3K27ac",]
plot<-ggplot(HUVEC_strengthened_H3K27ac, aes(x=data, y=signal, fill=data)) + geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=2)+
  theme_classic()+ ggtitle("Stripes strengthened in IMR90 relative to HUVEC cells")+
  scale_x_discrete(name="", breaks=waiver(), labels=waiver(), limits=NULL)+
  scale_y_continuous(name="Chip-seq signal", breaks=waiver(), labels=scientific_10, limits=NULL)+
  theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size=10,family = "Helvetica"))+
  stat_compare_means(method = "wilcox.test",comparisons = list(c("HUVEC","IMR90")),method.args=list(alternative = "greater"),paired=TRUE,
                     aes(label =  paste0(..method.., 'p =', ..p.format..)), label.y = 150000)+
  theme(legend.position="top",axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=.5))
print(plot)
ggsave(filename=paste(folder,"HUVEC_H3K27ac.svg",sep=""),dpi = 300,width=5,height=5)
print(dev.off())

HUVEC_strengthened_CTCF=HUVEC_strengthened[HUVEC_strengthened$protien=="CTCF",]
plot<-ggplot(HUVEC_strengthened_CTCF, aes(x=data, y=signal, fill=data)) + geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=2)+
  theme_classic()+ ggtitle("Stripes strengthened in IMR90 relative to HUVEC cells")+
  scale_x_discrete(name="", breaks=waiver(), labels=waiver(), limits=NULL)+
  scale_y_continuous(name="Chip-seq signal", breaks=waiver(), labels=scientific_10, limits=NULL)+
  theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size=10,family = "Helvetica"))+
  stat_compare_means(method = "wilcox.test",comparisons = list(c("HUVEC","IMR90")),method.args=list(alternative = "greater"),paired=TRUE,
                     aes(label =  paste0(..method.., 'p =', ..p.format..)), label.y = 150000)+
  theme(legend.position="top",axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=.5))
print(plot)
ggsave(filename=paste(folder,"HUVEC_CTCF.svg",sep=""),dpi = 300,width=5,height=5)
print(dev.off())



IMR90_strengthened=c()
for(i in 1:dim(IMR90)[1]){
  for(j in 1:length(IMR90[i,])){
    IMR90_strengthened=rbind(IMR90_strengthened,c(IMR90[i,][j][[1]],c(strsplit(names(IMR90[i,][j]),"_")[[1]]),i))
  }
}
colnames(IMR90_strengthened)<-c("signal","protien","data","id")
IMR90_strengthened=as.data.frame(IMR90_strengthened)
IMR90_strengthened$signal=as.numeric(IMR90_strengthened$signal)

IMR90_strengthened_Pol2=IMR90_strengthened[IMR90_strengthened$protien=="Pol2",]
plot<-ggplot(IMR90_strengthened_Pol2, aes(x=data, y=signal, fill=data)) + geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=2)+
  theme_classic()+ ggtitle("Stripes strengthened in HUVEC relative to IMR90 cells")+
  scale_x_discrete(name="", breaks=waiver(), labels=waiver(), limits=NULL)+
  scale_y_continuous(name="Chip-seq signal", breaks=waiver(), labels=scientific_10, limits=NULL)+
  theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size=10,family = "Helvetica"))+
  stat_compare_means(method = "wilcox.test",comparisons = list(c("IMR90","HUVEC")),method.args=list(alternative = "greater"),paired=TRUE,
                     aes(label =  paste0(..method.., 'p =', ..p.format..)), label.y = 750)+
  theme(legend.position="top",axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=.5))
print(plot)
ggsave(filename=paste(folder,"IMR90_Pol2.svg",sep=""),dpi = 300,width=5,height=5)
print(dev.off())

IMR90_strengthened_H3K27ac=IMR90_strengthened[IMR90_strengthened$protien=="H3K27ac",]
plot<-ggplot(IMR90_strengthened_H3K27ac, aes(x=data, y=signal, fill=data)) + geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=2)+
  theme_classic()+ ggtitle("Stripes strengthened in HUVEC relative to IMR90 cells")+
  scale_x_discrete(name="", breaks=waiver(), labels=waiver(), limits=NULL)+
  scale_y_continuous(name="Chip-seq signal", breaks=waiver(), labels=scientific_10, limits=NULL)+
  theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size=10,family = "Helvetica"))+
  stat_compare_means(method = "wilcox.test",comparisons = list(c("IMR90","HUVEC")),method.args=list(alternative = "greater"),paired=TRUE,
                     aes(label =  paste0(..method.., 'p =', ..p.format..)), label.y = 150000)+
  theme(legend.position="top",axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=.5))
print(plot)
ggsave(filename=paste(folder,"IMR90_H3K27ac.svg",sep=""),dpi = 300,width=5,height=5)
print(dev.off())

IMR90_strengthened_CTCF=IMR90_strengthened[IMR90_strengthened$protien=="CTCF",]
plot<-ggplot(IMR90_strengthened_CTCF, aes(x=data, y=signal, fill=data)) + geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=2)+
  theme_classic()+ ggtitle("Stripes strengthened in HUVEC relative to IMR90 cells")+
  scale_x_discrete(name="", breaks=waiver(), labels=waiver(), limits=NULL)+
  scale_y_continuous(name="Chip-seq signal", breaks=waiver(), labels=scientific_10, limits=NULL)+
  theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size=10,family = "Helvetica"))+
  stat_compare_means(method = "wilcox.test",comparisons = list(c("IMR90","HUVEC")),method.args=list(alternative = "greater"),paired=TRUE,
                     aes(label =  paste0(..method.., 'p =', ..p.format..)), label.y = 150000)+
  theme(legend.position="top",axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=.5))
print(plot)
ggsave(filename=paste(folder,"IMR90_CTCF.svg",sep=""),dpi = 300,width=5,height=5)
print(dev.off())