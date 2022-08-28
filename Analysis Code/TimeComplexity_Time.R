homeDir='~/Revision/1_6'
dataDir="/Revision/test_updated"

library(ggplot2)
library(doParallel)
library(foreach)



size_list=c(1000000,100000,50000,25000,10000)[1:5]
# type_list=c("GM12878_primary","GM12878_replicate")
# type_list=c("IMR90","HUVEC")
# type_list=c("IMR90","GM12878_primary")
# type_list=c("IMR90","GM12878_replicate")
id_list=c(paste("chr",1:22,sep=""),"chrX")[1:23]
dataList=list(c("GM12878_primary","GM12878_replicate"),c("IMR90","HUVEC"),c("IMR90","GM12878_primary"),c("IMR90","GM12878_replicate"))
ddd=c(1:4)[1:4]
w=300


timeCount=matrix(0,ncol=9,nrow=0)
fname=c("cptMean","score","rollapply","loess","argmax","diffStripe","localPeak","matchPeaks","edging",
        "bcpMean","stripeCalling","getStripe","stripeLocationAndLength","PvalScore","stripeDiffComplete")[c(10:15)]
cname=c("timeSum","timeMean","timeFrequency","functions","chromosome","bin","res","window","data")
colnames(timeCount)=cname
timeCount=data.frame(timeCount)

size_list=c(1000000,100000,50000,25000,10000)[1:5]
# type_list=c("GM12878_primary","GM12878_replicate")
type_list=c("IMR90_HUVEC","GM12878_primary_GM12878_replicate","IMR90_GM12878_replicate","IMR90_GM12878_primary")
id_list=c(paste("chr",1:22,sep=""),"chrX")[1:23]
w=300

for(binSize in size_list){
  for(type in type_list){
    for(id in id_list){
      if(binSize<1000000){
        res=paste(binSize/1000,"kb",sep="")
      }else{
        res=paste(binSize/1000000,"mb",sep="")
      }
      if(binSize>100000){
        w=300
      }else if(binSize>10000){
        w=300
      }else{
        w=300
      }
      if(id=="chr1"){
        timeSum=c(rep(0,length(fname)))
        names(timeSum)=fname
        timeMean=c(rep(0,length(fname)))
        names(timeMean)=fname
        timeFreq=c(rep(0,length(fname)))
        names(timeFreq)=fname
      }
      for(fileN in fname){
        folderN=paste(homeDir,"/files/Time/",type,"_",id,"_",res,"_","obs_VC1/",sep="")
        if(file.exists(paste(folderN,fileN,".txt",sep=""))){
          print(paste(type,id,res,w,fileN,"VC",timeDate::timeDate(),sep="_"))
          timeC=read.table(paste(folderN,fileN,".txt",sep=""))
          if(fileN=="bcpMean"){
            realTime=as.numeric(timeC$x) # just to remember Edging(stripeCalling) time from getStripe
            lastStripeCallingTime1=realTime
          }else if(fileN=="stripeCalling"){
            realTime=as.numeric(timeC$x)
            realTime=realTime -(sum(as.numeric(lastStripeCallingTime1))/length(realTime)) # just to remember Edging(stripeCalling) time from getStripe
            lastStripeCallingTime2=realTime
          }else if(fileN=="getStripe"){
            realTime=as.numeric(timeC$x)
            realTime=rep(realTime,4)
            realTime=realTime -(sum(as.numeric(lastStripeCallingTime2))/length(realTime)) # just for subtracting Edging(stripeCalling) time from getStripe
            print(paste(sum(lastStripeCallingTime1),sum(lastStripeCallingTime2),sum(realTime),sep="_"))
          }else if(fileN=="stripeLocationAndLength"){
            realTime=as.numeric(timeC$x)
            realTime=rep(realTime,2)
          }else{
            realTime=as.numeric(timeC$x)
          }
          temp1=matrix(c(sum(as.numeric(realTime)),mean(as.numeric(realTime)),length(realTime),fileN,id,binSize,res,w,type),ncol=length(cname))
          colnames(temp1)=cname
          if(sum(is.na(temp1))==0){
            print("done")
            timeCount=rbind(timeCount,temp1)
            timeSum[fileN]=timeSum[fileN]+sum(as.numeric(realTime))
            timeMean[fileN]=timeMean[fileN]+mean(as.numeric(realTime))
            timeFreq[fileN]=timeFreq[fileN]+length(as.numeric(realTime))
          }
        }
        if(id=="chrX"){
          print("Done")
          temp2=c(timeSum[fileN],timeMean[fileN]/length(id_list),timeFreq[fileN],fileN,"genome",binSize,res,w,type)
          names(temp2)=cname
          timeCount=rbind(timeCount,temp2)
        }
      }
    }
  }
}

timeCount1=timeCount


timeCount1[timeCount1$data=="GM12878_primary_GM12878_replicate" & timeCount1$chromosome=="genome",]
timeCount1[timeCount1$data=="IMR90_HUVEC" & timeCount1$chromosome=="genome",]

timeCount1=timeCount1[timeCount1$chromosome=="genome",]
timeCount1$data[timeCount1$data=="GM12878_primary_GM12878_replicate"]<-"GM12878_primary_replicate"
timeCount1=timeCount1[timeCount1$functions!="bcpMean",]
timeCount1$functions[timeCount1$functions=="stripeCalling"]<-"1. Edging"
timeCount1$functions[timeCount1$functions=="getStripe"]<-"2. stripe-detect"
timeCount1$functions[timeCount1$functions=="stripeLocationAndLength"]<-"3. stripe-length"
timeCount1$functions[timeCount1$functions=="PvalScore"]<-"4. Chi-square"
timeCount1$functions[timeCount1$functions=="stripeDiffComplete"]<-"Total 1-4"

# timeCount1$timeSum=round(as.numeric(timeCount1$timeSum)/60,2)
timeCount1$timeSum=round(log2((as.numeric(timeCount1$timeSum)/60)+1),2)

timeCount1$timeMean=round(log2((as.numeric(timeCount1$timeMean)/60)+1),2)



plot=ggplot(timeCount1, aes(x=res, y=timeSum, group=functions)) + 
  geom_point(aes(color=functions),size = 1)+
  geom_line(aes(color=functions),size = .2)+  
  theme_classic()+  ggtitle(paste("",sep=""))+ # theme_linedraw()+
  scale_x_discrete(name="Resolution", breaks=waiver(), labels=waiver(), limits=c("1mb","100kb","50kb","25kb","10kb"))+
  # scale_y_continuous(name="log time in minutes", breaks=c(2,4,6,8), labels=waiver(), limits=c(0,10))+
  scale_y_continuous(name="Time in minutes", breaks=c(2,4,6,8,10,12), labels=c(4,16,64,256,1024,4096), limits=c(0,12))+
  theme(text = element_text(size=10,family = "Helvetica"))+
  theme(legend.position="right",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))# +coord_flip()
plot+facet_grid(~data)
ggsave(filename=paste(homeDir,"/Time_Latest.pdf",sep=""),dpi = 300,width=5,height=5)
dev.off()



