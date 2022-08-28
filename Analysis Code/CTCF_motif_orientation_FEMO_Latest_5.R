library(ggplot2)
library(export)
library(Hmisc)
library(doParallel)
library(foreach)
home="~/Revision"

tt="VC1"
depth_list=c("")
cl=length(depth_list)
registerDoParallel(cl)
# for(depth in depth_list){
foreach::foreach(depth=depth_list,.combine=cbind)%dopar%{
  for (pval in c(.05)){
    folder=paste(home,"/2_4_1/data_Latest_5",depth,"_",tt,"_",pval,sep="")
    folder1=paste(home,"/2_4_1/data",depth,"_",tt,sep="")
    if(!dir.exists(folder)){
      dir.create(folder)
    }
    if(!dir.exists(folder1)){
      dir.create(folder1)
    }
    # direction left means, horihental stripe in upper right traingle
    # direction right means, horihental stripe in lower left traingle
    id_list=c(paste("chr",1:22,sep=""),"chrX")
    binSizeList=c(10000)
    for (binSize in binSizeList) {
      # binSize=1000000
      if(binSize<1000000){
        res=paste(binSize/1000,"kb",sep="")
      }else{
        res=paste(binSize/1000000,"mb",sep="")
      }
      
      
      type_list=c("IMR90","HUVEC")
      print(paste(res,pval,sep="_"))
      dataName=paste("A_5_stripe_",type_list[1],"_",type_list[2],"_",pval,sep="")
      if(file.exists(paste(folder,"/",dataName,"_",res,".rds",sep=""))){
        stripes1=readRDS(paste(folder,"/",dataName,"_",res,".rds",sep=""))
      }else{
        stripes1=list()
      }
      dataName1=paste("A_5_combine_",type_list[1],"_",type_list[2],sep="")
      if(file.exists(paste(folder,"/",dataName1,"_",res,".rds",sep=""))){
        combine=readRDS(paste(folder,"/",dataName1,"_",res,".rds",sep=""))
      }else{
        combine=list()
      }
      dataName2=paste("A_5_FinalResult_",type_list[1],"_",type_list[2],"_",pval,sep="")
      if(file.exists(paste(folder,"/",dataName2,"_",res,".rds",sep=""))){
        finalResult=readRDS(paste(folder,"/",dataName2,"_",res,".rds",sep=""))
      }else{
        finalResult=list()
      }
      
      stripes=c()
      flag=0
      for(id in id_list){
        print(id)
        files=paste(home,"/test_updated",depth,"/output_updated_",tt,"/",type_list[1],"_",type_list[2],depth,"_",
                    id,"_",res,"_obs_",tt,"/deduplicated_in_",type_list[1],"_not_",type_list[2],"_",id,"_stripes.txt",sep="")
        if(file.exists(files)){
          temp=read.table(files)
          if(flag==0){
            colname=c("in","out","res","comp",temp[1,])
            flag=1
          }
          stripes=rbind(stripes,cbind("in"=rep(type_list[1],dim(temp)[1]-1),"out"=rep(type_list[2],dim(temp)[1]-1),
                                      "res"=rep(res,dim(temp)[1]-1),"comp"=rep(paste(type_list[1],"_",type_list[2],sep=""),dim(temp)[1]-1),temp[-1,]))
        }
        files=paste(home,"/test_updated",depth,"/output_updated_",tt,"/",type_list[1],"_",type_list[2],depth,"_",
                    id,"_",res,"_obs_",tt,"/deduplicated_in_",type_list[2],"_not_",type_list[1],"_",id,"_stripes.txt",sep="")
        if(file.exists(files)){
          temp=read.table(files)
          stripes=rbind(stripes,cbind("in"=rep(type_list[2],dim(temp)[1]-1),"out"=rep(type_list[1],dim(temp)[1]-1),
                                      "res"=rep(res,dim(temp)[1]-1),"comp"=rep(paste(type_list[1],"_",type_list[2],sep=""),dim(temp)[1]-1),temp[-1,]))
        }
      }
      colnames(stripes)<-colname
      stripes<-na.omit(stripes)
      dim(stripes)
      
      for(id in id_list){
        dataL=paste("data_",type_list[1],"_",type_list[2],"_",id,sep="")
        if(file.exists(paste(folder1,"/",dataL,"_",res,".rds",sep=""))){
          dataList=readRDS(paste(folder1,"/",dataL,"_",res,".rds",sep=""))
        }else{
          dataList=list()
        }
        if(length(dataList)<2){
          print(id)
          dataList[[paste(type_list[1],"_",id,"_",res,"_obs_",tt,sep="")]]=as.matrix(readRDS(paste(home,"/test_updated",depth,"/files_",tt,"/",type_list[1],depth,"_",id,"_",res,"_obs_",tt,".rds",sep="")))
          dataList[[paste(type_list[2],"_",id,"_",res,"_obs_",tt,sep="")]]=as.matrix(readRDS(paste(home,"/test_updated",depth,"/files_",tt,"/",type_list[2],depth,"_",id,"_",res,"_obs_",tt,".rds",sep="")))
          saveRDS(dataList,paste(folder1,"/",dataL,"_",res,".rds",sep=""))
        }
        
        
        for(k in 1:dim(stripes)[1]){
          if(paste(stripes[k,"in"],"_",stripes[k,"chrom"],"_",res,"_obs_",tt,sep="") %in% names(dataList)){
            if(stripes[k,"in"]=="IMR90"){
              a=dataList[[paste(stripes[k,"in"],"_",stripes[k,"chrom"],"_",res,"_obs_",tt,sep="")]]
              ab=t(a)
              diag(ab)<-0
              a=a+ab
              b=dataList[[paste(stripes[k,"out"],"_",stripes[k,"chrom"],"_",res,"_obs_",tt,sep="")]]
              ab=t(b)
              diag(ab)<-0
              b=b+ab
            }else{
              a=dataList[[paste(stripes[k,"out"],"_",stripes[k,"chrom"],"_",res,"_obs_",tt,sep="")]]
              ab=t(a)
              diag(ab)<-0
              a=a+ab
              b=dataList[[paste(stripes[k,"in"],"_",stripes[k,"chrom"],"_",res,"_obs_",tt,sep="")]]
              ab=t(b)
              diag(ab)<-0
              b=b+ab
            }
            if(stripes[k,'direction']=="left"){
              if(as.numeric(stripes[k,"stripeLength"])<0){
                stripes[k,"stripeLength"]=0
              }
              fc=log2(mean(a[as.numeric(stripes[k,"leftEdge"]):as.numeric(stripes[k,"rightEdge"]),as.numeric(stripes[k,"leftEdge"]):as.numeric(stripes[k,"rightEdge"])])/
                        mean(b[as.numeric(stripes[k,"leftEdge"]):as.numeric(stripes[k,"rightEdge"]),as.numeric(stripes[k,"leftEdge"]):as.numeric(stripes[k,"rightEdge"])])) # saved in data with pval
              wilcoxPval=as.numeric(stripes[k,"p.chip"]) # Pval have been taken from StripeDiff method
              if(wilcoxPval<=pval){
                if(fc>=0){
                  sample1StripeType=paste(strsplit(as.character(stripes[k,"comp"]),"_")[[1]][1],'strengthened',sep="_")
                }else if(fc<=-0){
                  sample1StripeType=paste(strsplit(as.character(stripes[k,"comp"]),"_")[[1]][1],'weakened',sep="_")
                }
              }else{
                sample1StripeType='insignificant'
              }
              x=c("loc"=paste(stripes[k,"chrom"],":",as.numeric(stripes[k,"leftEdge"])*binSize,"-",binSize*as.numeric(stripes[k,"rightEdge"]),sep=""),
                  'firstLoc'=as.numeric(stripes[k,"leftEdge"])*binSize,'secondLoc'=binSize*as.numeric(stripes[k,"rightEdge"]),
                  "logFC"=fc,"wilcoxPval"=wilcoxPval,"sample1StripeType"=sample1StripeType,stripes[k,])
              if(length(stripes1)<1 || !x[["loc"]] %in% unlist(stripes1[,"loc"])){
                stripes1=rbind(stripes1,x)
                print(paste(dim(stripes)[1],"_",k,sep=""))
              }
            }else{
              if(as.numeric(stripes[k,"stripeLength"])<0){
                stripes[k,"stripeLength"]=0
              }
          
              fc=log2(mean(a[as.numeric(stripes[k,"leftEdge"]):as.numeric(stripes[k,"rightEdge"]),as.numeric(stripes[k,"leftEdge"]):as.numeric(stripes[k,"rightEdge"])])/
                        mean(b[as.numeric(stripes[k,"leftEdge"]):as.numeric(stripes[k,"rightEdge"]),as.numeric(stripes[k,"leftEdge"]):as.numeric(stripes[k,"rightEdge"])])) # saved in data with pval
              wilcoxPval=as.numeric(stripes[k,"p.chip"])
              if(wilcoxPval<=pval){
                if(fc>=0){
                  sample1StripeType=paste(strsplit(as.character(stripes[k,"comp"]),"_")[[1]][1],'strengthened',sep="_")
                }else if(fc<=-0){
                  sample1StripeType=paste(strsplit(as.character(stripes[k,"comp"]),"_")[[1]][1],'weakened',sep="_")
                }
              }else{
                sample1StripeType='insignificant'
              }
              x=c("loc"=paste(stripes[k,"chrom"],":",as.numeric(stripes[k,"leftEdge"])*binSize,"-",binSize*as.numeric(stripes[k,"rightEdge"]),sep=""),
                  'firstLoc'=as.numeric(stripes[k,"leftEdge"])*binSize,'secondLoc'=binSize*as.numeric(stripes[k,"rightEdge"]),
                  "logFC"=fc,"wilcoxPval"=wilcoxPval,"sample1StripeType"=sample1StripeType,stripes[k,])
              if(length(stripes1)<1 || !x[["loc"]] %in% unlist(stripes1[,"loc"])){
                stripes1=rbind(stripes1,x)
                print(paste(dim(stripes)[1],"_",k,sep=""))
              }
            }
          }
        }
        saveRDS(stripes1,paste(folder,"/",dataName,"_",res,".rds",sep=""))
      }
      
      
      stripes1=readRDS(paste(folder,"/",dataName,"_",res,".rds",sep=""))
      stripes1=data.frame(stripes1)
      table(unlist(stripes1$sample1StripeType))
      dim(stripes1)
      ctcf=read.table(paste(folder,"/../",'fimo.tsv',sep=''))
      dim(ctcf)
      colnames(ctcf)=ctcf[1,]
      ctcf=ctcf[-1,]
      ######################################3
      ctcf=ctcf[ctcf$`q-value`<=pval,]
      dim(ctcf)
      if((length(combine)+1)<=dim(stripes1)[1]){
        for(i in (length(combine)+1):dim(stripes1)[1]){
          if(!is.na(stripes1[i,c('loc')][[1]])){
            combine[[stripes1[i,c('loc')][[1]]]]<-as.list(as.matrix(stripes1[i,])[1,])
            combine[[stripes1[i,c('loc')][[1]]]][["count"]]=0
            combine[[stripes1[i,c('loc')][[1]]]][["count+"]]=0
            combine[[stripes1[i,c('loc')][[1]]]][["count-"]]=0
            combine[[stripes1[i,c('loc')][[1]]]][["randomCount+"]]=0
            combine[[stripes1[i,c('loc')][[1]]]][["randomCount-"]]=0
            combine[[stripes1[i,c('loc')][[1]]]][["leftcount+"]]=0
            combine[[stripes1[i,c('loc')][[1]]]][["leftcount-"]]=0
            combine[[stripes1[i,c('loc')][[1]]]][["rightcount+"]]=0
            combine[[stripes1[i,c('loc')][[1]]]][["rightcount-"]]=0
            combine[[stripes1[i,c('loc')][[1]]]][["leftcount"]]=0
            combine[[stripes1[i,c('loc')][[1]]]][["rightcount"]]=0
            combine[[stripes1[i,c('loc')][[1]]]][["strand"]]=c()
            combine[[stripes1[i,c('loc')][[1]]]][["strand+"]]=c()
            combine[[stripes1[i,c('loc')][[1]]]][["strand-"]]=c()
            combine[[stripes1[i,c('loc')][[1]]]][["randomStrand+"]]=c()
            combine[[stripes1[i,c('loc')][[1]]]][["randomStrand-"]]=c()
            combine[[stripes1[i,c('loc')][[1]]]][["leftstrand+"]]=c()
            combine[[stripes1[i,c('loc')][[1]]]][["leftstrand-"]]=c()
            combine[[stripes1[i,c('loc')][[1]]]][["rightstrand+"]]=c()
            combine[[stripes1[i,c('loc')][[1]]]][["rightstrand-"]]=c()
            combine[[stripes1[i,c('loc')][[1]]]][["leftstrand"]]=c()
            combine[[stripes1[i,c('loc')][[1]]]][["rightstrand"]]=c()
            combine[[stripes1[i,c('loc')][[1]]]][["start"]]=c()
            combine[[stripes1[i,c('loc')][[1]]]][["stop"]]=c()
            for(j in 1:dim(ctcf)[1]){
              print(paste(dim(stripes1)[1],i,j,sep="__"))
              stripeLocationDiff=as.numeric(stripes1[i,'secondLoc'][[1]])-as.numeric(stripes1[i,'firstLoc'][[1]])
              randomLocation=sample(1:(as.numeric(stripes1[i,'secondLoc'][[1]])*10),1)
              randomLocation1=randomLocation
              randomLocation2=randomLocation+stripeLocationDiff
              if(((as.numeric(ctcf[j,]$start) >= randomLocation1 & as.numeric(ctcf[j,]$start) <= randomLocation2) |
                  (as.numeric(ctcf[j,]$stop) >= randomLocation1 & as.numeric(ctcf[j,]$stop) <= randomLocation2)) &
                 ctcf[j,]$sequence_name==stripes1[i,"chrom"][[1]]){ # fast checking condition
                if(ctcf[j,]$strand=="+"){
                  combine[[stripes1[i,c('loc')][[1]]]][["randomCount+"]]=combine[[stripes1[i,c('loc')][[1]]]][["randomCount+"]]+1
                  combine[[stripes1[i,c('loc')][[1]]]][["randomStrand+"]]=c(combine[[stripes1[i,c('loc')][[1]]]][["randomStrand+"]],ctcf[j,]$strand)
                }
                if(ctcf[j,]$strand=="-"){
                  combine[[stripes1[i,c('loc')][[1]]]][["randomCount-"]]=combine[[stripes1[i,c('loc')][[1]]]][["randomCount-"]]+1
                  combine[[stripes1[i,c('loc')][[1]]]][["randomStrand-"]]=c(combine[[stripes1[i,c('loc')][[1]]]][["randomStrand-"]],ctcf[j,]$strand)
                }
              }
              
              if(((as.numeric(ctcf[j,]$start) >= as.numeric(stripes1[i,'firstLoc'][[1]]) & as.numeric(ctcf[j,]$start) <= as.numeric(stripes1[i,'secondLoc'][[1]])) |
                  (as.numeric(ctcf[j,]$stop) >= as.numeric(stripes1[i,'firstLoc'][[1]]) & as.numeric(ctcf[j,]$stop) <= as.numeric(stripes1[i,'secondLoc'][[1]]))) &
                 ctcf[j,]$sequence_name==stripes1[i,"chrom"][[1]]){ # fast checking condition
                combine[[stripes1[i,c('loc')][[1]]]][["count"]]=combine[[stripes1[i,c('loc')][[1]]]][["count"]]+1
                combine[[stripes1[i,c('loc')][[1]]]][["strand"]]=c(combine[[stripes1[i,c('loc')][[1]]]][["strand"]],ctcf[j,]$strand)
                if(ctcf[j,]$strand=="+"){
                  combine[[stripes1[i,c('loc')][[1]]]][["count+"]]=combine[[stripes1[i,c('loc')][[1]]]][["count+"]]+1
                  combine[[stripes1[i,c('loc')][[1]]]][["strand+"]]=c(combine[[stripes1[i,c('loc')][[1]]]][["strand+"]],ctcf[j,]$strand)
                }
                if(ctcf[j,]$strand=="-"){
                  combine[[stripes1[i,c('loc')][[1]]]][["count-"]]=combine[[stripes1[i,c('loc')][[1]]]][["count-"]]+1
                  combine[[stripes1[i,c('loc')][[1]]]][["strand-"]]=c(combine[[stripes1[i,c('loc')][[1]]]][["strand-"]],ctcf[j,]$strand)
                }
                
                if(stripes1[i,"direction"][[1]]=="left"){
                  combine[[stripes1[i,c('loc')][[1]]]][["leftcount"]]=combine[[stripes1[i,c('loc')][[1]]]][["leftcount"]]+1
                  combine[[stripes1[i,c('loc')][[1]]]][["leftstrand"]]=c(combine[[stripes1[i,c('loc')][[1]]]][["leftstrand"]],ctcf[j,]$strand)
                  if(ctcf[j,]$strand=="+"){
                    combine[[stripes1[i,c('loc')][[1]]]][["leftcount+"]]=combine[[stripes1[i,c('loc')][[1]]]][["leftcount+"]]+1
                    combine[[stripes1[i,c('loc')][[1]]]][["leftstrand+"]]=c(combine[[stripes1[i,c('loc')][[1]]]][["leftstrand+"]],ctcf[j,]$strand)
                  }
                  if(ctcf[j,]$strand=="-"){
                    combine[[stripes1[i,c('loc')][[1]]]][["leftcount-"]]=combine[[stripes1[i,c('loc')][[1]]]][["leftcount-"]]+1
                    combine[[stripes1[i,c('loc')][[1]]]][["leftstrand-"]]=c(combine[[stripes1[i,c('loc')][[1]]]][["leftstrand-"]],ctcf[j,]$strand)
                  }
                }
                
                if(stripes1[i,"direction"][[1]]=="right"){
                  combine[[stripes1[i,c('loc')][[1]]]][["rightcount"]]=combine[[stripes1[i,c('loc')][[1]]]][["rightcount"]]+1
                  combine[[stripes1[i,c('loc')][[1]]]][["rightstrand"]]=c(combine[[stripes1[i,c('loc')][[1]]]][["rightstrand"]],ctcf[j,]$strand)
                  if(ctcf[j,]$strand=="+"){
                    combine[[stripes1[i,c('loc')][[1]]]][["rightcount+"]]=combine[[stripes1[i,c('loc')][[1]]]][["rightcount+"]]+1
                    combine[[stripes1[i,c('loc')][[1]]]][["rightstrand+"]]=c(combine[[stripes1[i,c('loc')][[1]]]][["rightstrand+"]],ctcf[j,]$strand)
                  }
                  if(ctcf[j,]$strand=="-"){
                    combine[[stripes1[i,c('loc')][[1]]]][["rightcount-"]]=combine[[stripes1[i,c('loc')][[1]]]][["rightcount-"]]+1
                    combine[[stripes1[i,c('loc')][[1]]]][["rightstrand-"]]=c(combine[[stripes1[i,c('loc')][[1]]]][["rightstrand-"]],ctcf[j,]$strand)
                  }
                }
                
                combine[[stripes1[i,c('loc')][[1]]]][["start"]]=c(combine[[stripes1[i,c('loc')][[1]]]][["start"]],as.numeric(ctcf[j,]$start))
                combine[[stripes1[i,c('loc')][[1]]]][["stop"]]=c(combine[[stripes1[i,c('loc')][[1]]]][["stop"]],as.numeric(ctcf[j,]$stop))
              }
            }
            if(length(combine[[stripes1[i,c('loc')][[1]]]][["strand"]])>0){
              combine[[stripes1[i,c('loc')][[1]]]][["motifPercentage1"]]=sum( (combine[[stripes1[i,c('loc')][[1]]]][["strand"]]=="-" & combine[[stripes1[i,c('loc')][[1]]]][["direction"]]=="left") |
                                                                                (combine[[stripes1[i,c('loc')][[1]]]][["strand"]]=="+" & combine[[stripes1[i,c('loc')][[1]]]][["direction"]]=="right") )/
                length(combine[[stripes1[i,c('loc')][[1]]]][["strand"]])*100
              
              combine[[stripes1[i,c('loc')][[1]]]][["motifPercentage2"]]=sum( (combine[[stripes1[i,c('loc')][[1]]]][["strand"]]=="+" & combine[[stripes1[i,c('loc')][[1]]]][["direction"]]=="left") |
                                                                                (combine[[stripes1[i,c('loc')][[1]]]][["strand"]]=="-" & combine[[stripes1[i,c('loc')][[1]]]][["direction"]]=="right") )/
                length(combine[[stripes1[i,c('loc')][[1]]]][["strand"]])*100
              
              combine[[stripes1[i,c('loc')][[1]]]][["motifPercentage3"]]=sum( (combine[[stripes1[i,c('loc')][[1]]]][["strand"]]=="-" & combine[[stripes1[i,c('loc')][[1]]]][["direction"]]=="left") |
                                                                                (combine[[stripes1[i,c('loc')][[1]]]][["strand"]]=="+" & combine[[stripes1[i,c('loc')][[1]]]][["direction"]]=="right") )/
                sum( (combine[[stripes1[i,c('loc')][[1]]]][["strand"]]=="-" & combine[[stripes1[i,c('loc')][[1]]]][["direction"]]=="left") |
                       (combine[[stripes1[i,c('loc')][[1]]]][["strand"]]=="+" & combine[[stripes1[i,c('loc')][[1]]]][["direction"]]=="right") |
                       (combine[[stripes1[i,c('loc')][[1]]]][["strand"]]=="+" & combine[[stripes1[i,c('loc')][[1]]]][["direction"]]=="left") |
                       (combine[[stripes1[i,c('loc')][[1]]]][["strand"]]=="-" & combine[[stripes1[i,c('loc')][[1]]]][["direction"]]=="right"))*100
              
              combine[[stripes1[i,c('loc')][[1]]]][["motifPercentage4"]]=sum( (combine[[stripes1[i,c('loc')][[1]]]][["strand"]]=="+" & combine[[stripes1[i,c('loc')][[1]]]][["direction"]]=="left") |
                                                                                (combine[[stripes1[i,c('loc')][[1]]]][["strand"]]=="-" & combine[[stripes1[i,c('loc')][[1]]]][["direction"]]=="right") )/
                sum( (combine[[stripes1[i,c('loc')][[1]]]][["strand"]]=="-" & combine[[stripes1[i,c('loc')][[1]]]][["direction"]]=="left") |
                       (combine[[stripes1[i,c('loc')][[1]]]][["strand"]]=="+" & combine[[stripes1[i,c('loc')][[1]]]][["direction"]]=="right") |
                       (combine[[stripes1[i,c('loc')][[1]]]][["strand"]]=="+" & combine[[stripes1[i,c('loc')][[1]]]][["direction"]]=="left") |
                       (combine[[stripes1[i,c('loc')][[1]]]][["strand"]]=="-" & combine[[stripes1[i,c('loc')][[1]]]][["direction"]]=="right"))*100
              
              combine[[stripes1[i,c('loc')][[1]]]][["motifPercentage5"]]=sum(combine[[stripes1[i,c('loc')][[1]]]][["strand"]]=="+")/
                length(combine[[stripes1[i,c('loc')][[1]]]][["strand"]])*100
              
              combine[[stripes1[i,c('loc')][[1]]]][["motifPercentage6"]]=sum(combine[[stripes1[i,c('loc')][[1]]]][["strand"]]=="-")/
                length(combine[[stripes1[i,c('loc')][[1]]]][["strand"]])*100
            }else{
              combine[[stripes1[i,c('loc')][[1]]]][["motifPercentage1"]]=0
              combine[[stripes1[i,c('loc')][[1]]]][["motifPercentage2"]]=0
              combine[[stripes1[i,c('loc')][[1]]]][["motifPercentage3"]]=0
              combine[[stripes1[i,c('loc')][[1]]]][["motifPercentage4"]]=0
              combine[[stripes1[i,c('loc')][[1]]]][["motifPercentage5"]]=0
              combine[[stripes1[i,c('loc')][[1]]]][["motifPercentage6"]]=0
            }
            saveRDS(combine,paste(folder,"/",dataName1,"_",res,".rds",sep=""))
          }
        }
      }
      
      # saveRDS(combine,paste(folder,"/",dataName1,"_",res,".rds",sep=""))
      combine=readRDS(paste(folder,"/",dataName1,"_",res,".rds",sep=""))
      
    }
  }
}



#####################combine plot
library(ggplot2)
library(ggpubr)
home="~/Desktop/stripeDiff-master/stripeDiff_Krishan/Revision"
# home="/home/ch228301/Public/Krishan/Krishan_Gupta/stripeDiff_master/stripeDiff_Krishan/Revision"
type_list=c("IMR90","HUVEC")
dataName1=paste("A_5_combine_",type_list[1],"_",type_list[2],sep="")
tt="VC1"
Allsets=c()
Allsets1=c()
ctcfBindingMotif=c()
ctcfBindingMotif1=c()
depth=""
res="10kb"
pval=0.05
id_min=1
id_1=1

folder=paste(home,"/2_4_1/data_Latest_5",depth,"_",tt,"_",pval,sep="")
binSize=10000
id_list=c(paste("chr",1:22,sep=""),"chrX")
combine=readRDS(paste(folder,"/",dataName1,"_",res,".rds",sep=""))
id_max=id_min+length(combine)-1
print(paste("first",depth,res,id_min,id_max,sep="_"))
id=id_min:id_max
id_min=id_max+1
combine1<-data.frame("StripeType"=unlist(lapply(combine,function(x) x$sample1StripeType)),
                     "motif1"=unlist(lapply(combine,function(x) x$motifPercentage1)),
                     "motif2"=unlist(lapply(combine,function(x) x$motifPercentage2)),
                     "motif3"=unlist(lapply(combine,function(x) x$motifPercentage3)),
                     "motif4"=unlist(lapply(combine,function(x) x$motifPercentage4)),
                     "motif5"=unlist(lapply(combine,function(x) x$motifPercentage5)),
                     "motif6"=unlist(lapply(combine,function(x) x$motifPercentage6)),
                     "both"=unlist(lapply(combine,function(x) x$count)),
                     "positive"=unlist(lapply(combine,function(x) x$`count+`)),
                     "negative"=unlist(lapply(combine,function(x) x$`count-`)),
                     "randomPositive"=unlist(lapply(combine,function(x) x$`randomCount+`)),
                     "randomNegative"=unlist(lapply(combine,function(x) x$`randomCount-`)),
                     "leftpositive"=unlist(lapply(combine,function(x) x$`leftcount+`)),
                     "leftnegative"=unlist(lapply(combine,function(x) x$`leftcount-`)),
                     "rightpositive"=unlist(lapply(combine,function(x) x$`rightcount+`)),
                     "rightnegative"=unlist(lapply(combine,function(x) x$`rightcount-`)),
                     "left"=unlist(lapply(combine,function(x) x$`leftcount`)),
                     "right"=unlist(lapply(combine,function(x) x$`rightcount`)),
                     "maxStripe"=unlist(lapply(combine,function(x) max(x$`count+`,x$`count-`))),
                     "minStripe"=unlist(lapply(combine,function(x) min(x$`count+`,x$`count-`))),
                     "maxRandom"=unlist(lapply(combine,function(x) max(x$`randomCount+`,x$`randomCount-`))),
                     "minRandom"=unlist(lapply(combine,function(x) min(x$`randomCount+`,x$`randomCount-`))),"ID"=id)

for(kkkk in 1:dim(combine1)[1]){
  print(paste(depth,res,kkkk,combine1$ID[kkkk],sep="_"))
  ctcfBindingMotif1=rbind(ctcfBindingMotif1,cbind("count"=combine1$negative[kkkk],
                                                  "strand"="-","ID"=combine1$ID[kkkk],
                                                  "StripeType"=combine1$StripeType[kkkk],"res"=res,"depth"=depth))
  ctcfBindingMotif1=rbind(ctcfBindingMotif1,cbind("count"=combine1$positive[kkkk],
                                                  "strand"="+","ID"=combine1$ID[kkkk],
                                                  "StripeType"=combine1$StripeType[kkkk],"res"=res,"depth"=depth))
  ctcfBindingMotif1=rbind(ctcfBindingMotif1,cbind("count"=combine1$randomNegative[kkkk],
                                                  "strand"="r-","ID"=combine1$ID[kkkk],
                                                  "StripeType"=combine1$StripeType[kkkk],"res"=res,"depth"=depth))
  ctcfBindingMotif1=rbind(ctcfBindingMotif1,cbind("count"=combine1$randomPositive[kkkk],
                                                  "strand"="r+","ID"=combine1$ID[kkkk],
                                                  "StripeType"=combine1$StripeType[kkkk],"res"=res,"depth"=depth))
  ctcfBindingMotif1=rbind(ctcfBindingMotif1,cbind("count"=combine1$leftnegative[kkkk],
                                                  "strand"="left-","ID"=combine1$ID[kkkk],
                                                  "StripeType"=combine1$StripeType[kkkk],"res"=res,"depth"=depth))
  ctcfBindingMotif1=rbind(ctcfBindingMotif1,cbind("count"=combine1$leftpositive[kkkk],
                                                  "strand"="left+","ID"=combine1$ID[kkkk],
                                                  "StripeType"=combine1$StripeType[kkkk],"res"=res,"depth"=depth))
  ctcfBindingMotif1=rbind(ctcfBindingMotif1,cbind("count"=combine1$rightnegative[kkkk],
                                                  "strand"="right-","ID"=combine1$ID[kkkk],
                                                  "StripeType"=combine1$StripeType[kkkk],"res"=res,"depth"=depth))
  ctcfBindingMotif1=rbind(ctcfBindingMotif1,cbind("count"=combine1$rightpositive[kkkk],
                                                  "strand"="right+","ID"=combine1$ID[kkkk],
                                                  "StripeType"=combine1$StripeType[kkkk],"res"=res,"depth"=depth))
  ctcfBindingMotif1=rbind(ctcfBindingMotif1,cbind("count"=combine1$left[kkkk],
                                                  "strand"="left","ID"=combine1$ID[kkkk],
                                                  "StripeType"=combine1$StripeType[kkkk],"res"=res,"depth"=depth))
  ctcfBindingMotif1=rbind(ctcfBindingMotif1,cbind("count"=combine1$right[kkkk],
                                                  "strand"="right","ID"=combine1$ID[kkkk],
                                                  "StripeType"=combine1$StripeType[kkkk],"res"=res,"depth"=depth))
  ctcfBindingMotif1=rbind(ctcfBindingMotif1,cbind("count"=combine1$both[kkkk],
                                                  "strand"="both","ID"=combine1$ID[kkkk],
                                                  "StripeType"=combine1$StripeType[kkkk],"res"=res,"depth"=depth))
  ctcfBindingMotif1=rbind(ctcfBindingMotif1,cbind("count"=sum(combine1$randomPositive[kkkk]+combine1$randomNegative[kkkk]),
                                                  "strand"="bothR","ID"=combine1$ID[kkkk],
                                                  "StripeType"=combine1$StripeType[kkkk],"res"=res,"depth"=depth))
}
    

# saveRDS(ctcfBindingMotif1,paste(home,"/2_4_1/ctcfBindingMotif1_5.rds",sep=""))

ctcfBindingMotif1=readRDS(paste(home,"/2_4_1/ctcfBindingMotif1_5.rds",sep=""))




#####################################################################
ctcfBindingMotif1=data.frame(ctcfBindingMotif1)
ctcfBindingMotif1=ctcfBindingMotif1[which(sort(ctcfBindingMotif1$ID) %in% ctcfBindingMotif1$ID),]
ctcfBindingMotif1$count=as.numeric(ctcfBindingMotif1$count)
ctcfBindingMotif1$strand=factor(ctcfBindingMotif1$strand,levels = c("both","left","right","-","+","left-","right-","left+","right+","r+","r-","bothR"))

ss=c("","left+","left-","right-","right+")
ctcfBindingMotif2=ctcfBindingMotif1[ctcfBindingMotif1$strand %in% ss,]
ctcfBindingMotif2$strand=factor(ctcfBindingMotif2$strand,levels = ss)


ctcfBindingMotif2=ctcfBindingMotif2[ctcfBindingMotif2$count>=0,]



standard_error <- function(x) sd(x) / sqrt(length(x))
ctcfBindingMotif3=ctcfBindingMotif2[ctcfBindingMotif2$depth %in% c(""),]
ctcfBindingMotif3=ctcfBindingMotif3[ctcfBindingMotif3$res %in% c("10kb"),]
ctcfBindingMotif3=rbind(cbind(ctcfBindingMotif3[ctcfBindingMotif3$strand %in% c("left+","right-"),],
                              "type"="stripe",
                              "stripeStrand"="left+, right-",
                              "orientation"="orientation",
                              "mean"=mean(ctcfBindingMotif3[ctcfBindingMotif3$strand %in% c("left+","right-"),"count"]),
                              "stderror"=standard_error(ctcfBindingMotif3[ctcfBindingMotif3$strand %in% c("left+","right-"),"count"]),
                              "Pval"=wilcox.test(ctcfBindingMotif3[ctcfBindingMotif3$strand %in% c("left+","right-"),"count"],
                                                 ctcfBindingMotif3[ctcfBindingMotif3$strand %in% c("left-","right+"),"count"],paired = TRUE,alternative = "greater")$p.value,
                              "ID_2"=c(ctcfBindingMotif3$ID[ctcfBindingMotif3$strand=="left+"],
                                       as.numeric(max(ctcfBindingMotif3$ID[ctcfBindingMotif3$strand=="left+"]))+
                                         as.numeric(ctcfBindingMotif3$ID[ctcfBindingMotif3$strand=="right-"]))),
                        cbind(ctcfBindingMotif3[ctcfBindingMotif3$strand %in% c("left-","right+"),],
                              "type"="stripe",
                              "stripeStrand"="left-, right+",
                              "orientation"="reverse",
                              "mean"=mean(ctcfBindingMotif3[ctcfBindingMotif3$strand %in% c("left-","right+"),"count"]),
                              "stderror"=standard_error(ctcfBindingMotif3[ctcfBindingMotif3$strand %in% c("left-","right+"),"count"]),
                              "Pval"=wilcox.test(ctcfBindingMotif3[ctcfBindingMotif3$strand %in% c("left+","right-"),"count"],
                                                 ctcfBindingMotif3[ctcfBindingMotif3$strand %in% c("left-","right+"),"count"],paired = TRUE,alternative = "greater")$p.value,
                              "ID_2"=c(ctcfBindingMotif3$ID[ctcfBindingMotif3$strand=="left-"],
                                       as.numeric(max(ctcfBindingMotif3$ID[ctcfBindingMotif3$strand=="left-"]))+
                                         as.numeric(ctcfBindingMotif3$ID[ctcfBindingMotif3$strand=="right+"]))))

ctcfBindingMotif3$ID_2=as.numeric(ctcfBindingMotif3$ID_2)
ctcfBindingMotifFinal1=ctcfBindingMotif3
######################################################################################################### differntial Stripe
standard_error <- function(x) sd(x) / sqrt(length(x))
ctcfBindingMotif3=ctcfBindingMotif2[ctcfBindingMotif2$depth %in% c(""),]
ctcfBindingMotif3=ctcfBindingMotif3[ctcfBindingMotif3$res %in% c("10kb"),]
ctcfBindingMotif3=ctcfBindingMotif3[!ctcfBindingMotif3$StripeType %in% c("insignificant"),]
ctcfBindingMotif3=rbind(cbind(ctcfBindingMotif3[ctcfBindingMotif3$strand %in% c("left+","right-"),],
                              "type"="diffStripe",
                              "stripeStrand"="left+, right-",
                              "orientation"="orientation",
                              "mean"=mean(ctcfBindingMotif3[ctcfBindingMotif3$strand %in% c("left+","right-"),"count"]),
                              "stderror"=standard_error(ctcfBindingMotif3[ctcfBindingMotif3$strand %in% c("left+","right-"),"count"]),
                              "Pval"=wilcox.test(ctcfBindingMotif3[ctcfBindingMotif3$strand %in% c("left+","right-"),"count"],
                                                 ctcfBindingMotif3[ctcfBindingMotif3$strand %in% c("left-","right+"),"count"],paired = TRUE,alternative = "greater")$p.value,
                              "ID_2"=c(ctcfBindingMotif3$ID[ctcfBindingMotif3$strand=="left+"],
                                       as.numeric(max(ctcfBindingMotif3$ID[ctcfBindingMotif3$strand=="left+"]))+
                                         as.numeric(ctcfBindingMotif3$ID[ctcfBindingMotif3$strand=="right-"]))),
                        cbind(ctcfBindingMotif3[ctcfBindingMotif3$strand %in% c("left-","right+"),],
                              "type"="diffStripe",
                              "stripeStrand"="left-, right+",
                              "orientation"="reverse",
                              "mean"=mean(ctcfBindingMotif3[ctcfBindingMotif3$strand %in% c("left-","right+"),"count"]),
                              "stderror"=standard_error(ctcfBindingMotif3[ctcfBindingMotif3$strand %in% c("left-","right+"),"count"]),
                              "Pval"=wilcox.test(ctcfBindingMotif3[ctcfBindingMotif3$strand %in% c("left+","right-"),"count"],
                                                 ctcfBindingMotif3[ctcfBindingMotif3$strand %in% c("left-","right+"),"count"],paired = TRUE,alternative = "greater")$p.value,
                              "ID_2"=c(ctcfBindingMotif3$ID[ctcfBindingMotif3$strand=="left-"],
                                       as.numeric(max(ctcfBindingMotif3$ID[ctcfBindingMotif3$strand=="left-"]))+
                                         as.numeric(ctcfBindingMotif3$ID[ctcfBindingMotif3$strand=="right+"]))))

ctcfBindingMotif3$ID_2=as.numeric(ctcfBindingMotif3$ID_2)
ctcfBindingMotifFinal2=ctcfBindingMotif3





library(export)
library(Hmisc)
ctcfBindingMotifFinal=rbind(ctcfBindingMotifFinal1,ctcfBindingMotifFinal2)
ctcfBindingMotifFinal=ctcfBindingMotifFinal[,c("res","depth","type","stripeStrand","orientation","mean","stderror","Pval")]
ctcfBindingMotifFinal=ctcfBindingMotifFinal[!duplicated(ctcfBindingMotifFinal),]

plot<-ggplot(ctcfBindingMotifFinal, aes(x=type, y=mean,fill=orientation)) +  
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar( aes(x=type, ymin=mean-stderror, ymax=mean+stderror), width=0.2, colour="black", alpha=0.9, size=1,position=position_dodge(.9)) +
  theme_classic()+ ggtitle("CTCF binding motifs orientation at 10kb resolution (IMR90 vs HUVEC)")+ 
  scale_x_discrete(name="", breaks=c("stripe","diffStripe"), labels=waiver(), limits=c("stripe","diffStripe"))+
  scale_y_continuous(name="Orientation's mean and standard error", breaks=waiver(), labels=waiver(), c(0,20))+
  theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size=10,family = "Helvetica"))+
  annotate(geom="text",x=c(1,2),y=.3,label=paste("p = ",unique(ctcfBindingMotifFinal$Pval),sep=""),color="black",size=5)+
  theme(legend.position="top",axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=.5))
print(plot)
ggsave(filename=paste(home,"/2_4_1/anchorOrientationLatest_5",type_list[1],"_",type_list[2],"_",tt,"_finalFigure.pdf",sep=""),dpi = 300,width=5,height=5)
print(dev.off())






