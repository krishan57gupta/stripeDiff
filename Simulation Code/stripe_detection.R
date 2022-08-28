########## pipline for comparing straps for two samples #########

getStrap <- function(contactMap, upLimit, window, outputPath='NA', outputName1='NA', 
                     outputName2='NA', outputName3='NA', outputName4='NA', 
                     startPoint=0, endPoint=0, p=0.8){
  # contactMap: a matrix for Hi-C contact map
  # upLimit: a number of up limit of Hi-C contact
  # outputPath: the output file path 
  # outputName1: the file name for heatmap of differential contact map
  # outputName2: the file name for bcp plot for all peaks
  # outputName3: the file name for bcp plot for up peaks
  # outputName4: the file name for bcp plot for down peaks
  # startPoint: the start point for the sub-contact map
  # endPoint: the end point for the sub-contact map
  # p: p-value cutoff for peak calling
  # output: a list of strap, probability of up peak, probability of down peak
  
  require(bcp)
  require(plotly)
  require(proxy)
  
  if(outputPath!='NA'){
    setwd(outputPath)
  }
  
  strapCalling <- function(contactMap, upLimit, outputPath, outputName1, outputName2, startPoint, endPoint, window){
    # input: Hi-C matrix
    # output: sum of differetial Hi-C matrix on each line
    
    sub_matrix <- function(map, w){
      for(i in 1:ncol(map)){
        if((i-w)>1){
          map[1:(i-w),i] = 0
        }
      }
      return(map)
    }
    
    contactMap = as.matrix(contactMap)
    contactMap[contactMap > upLimit] = upLimit
    contactMapUpperTri = contactMap
    contactMapUpperTri[lower.tri(contactMapUpperTri, diag = TRUE)] = 0
    contactMapUpperTri = sub_matrix(contactMapUpperTri, window)
    if(startPoint!=0&endPoint!=0){
      contactMapUpperTriMean = colMeans(t(diff(contactMapUpperTri)[startPoint:endPoint, startPoint:endPoint]))
    }
    else{
      contactMapUpperTriMean = colMeans(t(diff(contactMapUpperTri)))
    }
    
    bcpMean <- bcp(contactMapUpperTriMean,mcmc = 5000)
    
    if(outputName1!='NA'){
      if(startPoint!=0&endPoint!=0){
        s <- subplot(
          plot_ly(y = contactMapUpperTriMean, type = "bar"), 
          plot_ly(z = t(diff(contactMapUpperTri)[startPoint:endPoint, startPoint:endPoint]), type = "heatmap", 
                  zauto = F, zmin = -30, zmax = 20),
          nrows = 2, heights = c(0.5, 0.5), 
          shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = FALSE
        )}
      else{
        s <- subplot(
          plot_ly(y = contactMapUpperTriMean, type = "bar"), 
          plot_ly(z = t(diff(contactMapUpperTri)), type = "heatmap", zauto = F, zmin = -30, zmax = 20),
          nrows = 2, heights = c(0.5, 0.5), 
          shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = FALSE
        )
      }
      orca(s, outputName1)
    }
    
    if(outputName2!='NA'){
      pdf(outputName2, width = 4.5, height = 4)
      plot(bcpMean)
      dev.off()
    }
    return(bcpMean)
  }
  
  oneSideData <- function(bcpMean){
    # input: sum of diff Hi-C contact on each line
    # output: up peaks and down peaks
    rawData = bcpMean$data[,2]
    posteriorMeanLowProb = bcpMean$data[bcpMean$posterior.prob<0.1, 2]
    posteriorMeanLowProb = posteriorMeanLowProb[!is.na(posteriorMeanLowProb)]
    dens = density(posteriorMeanLowProb)
    sampleData= dens$x[dens$y>(max(dens$y)/2)]
    meanLowProb = mean(posteriorMeanLowProb)
    
    n1 = length(bcpMean$data[bcpMean$data[,2] < meanLowProb, 2])
    if(length(sampleData[sampleData<meanLowProb]) ==0){
      sampleDownSide = sample(meanLowProb, n1, replace = TRUE)
    }
    else{
      sampleDownSide = sample(sampleData[sampleData<meanLowProb], n1, replace = TRUE)
    }
    
    n2 = length(bcpMean$data[bcpMean$data[,2] >= meanLowProb, 2])
    if(length(sampleData[sampleData >= meanLowProb])==0){
      sampleUpSide = sample(meanLowProb, n2, replace = TRUE)
    }
    else{
      sampleUpSide = sample(sampleData[sampleData >= meanLowProb], n2, replace = TRUE)
    }
    
    upSideData = rawData
    upSideData[upSideData<meanLowProb] = sampleDownSide
    
    downSideData = rawData
    downSideData[downSideData>=meanLowProb] = sampleUpSide
    
    ########### use estimated TAD size to detect the down peaks ##############
    
    
    bcpUpSide = bcp(upSideData)
    bcpDownSide = bcp(downSideData)
    return(list(bcpUpSide, bcpDownSide))
  }
  
  matchPeaks <- function(upSide, downSide, p, map){
    upSideData <- upSide$data
    downSideData <- downSide$data
    upPeak = upSideData[upSide$posterior.prob>p, ]
    downPeak = downSideData[downSide$posterior.prob>p, ]
    
    if( !all(is.na(upPeak)) & !all(is.na(downPeak))){
      # cluster up-peaks and estimate the cutoff for distance between adjacent up peaks and down peaks
      boundaryDist = proxy::dist(upPeak[-nrow(upPeak),1], method = "euclidean")
      boundaryDist = as.matrix(boundaryDist)
      diagonal = boundaryDist[row(boundaryDist) == col(boundaryDist)+1]
      if(length(diagonal)==0){
        diagonal=1
      }
      if(length(unique(diagonal))==1){
        meanTadSize=diagonal}
      if(length(unique(diagonal))==2){
        cluster = kmeans(diagonal, 1)
        meanTadSize = max(cluster$centers)}
      if(length(unique(diagonal))>2){
        cluster = kmeans(diagonal, 2)
      meanTadSize = max(cluster$centers)}   # mean value for TAD size
      distance <- proxy::dist(data.frame(upPeak[,1]), data.frame(downPeak[,1]), method = 'euclidean')
      index = which(distance < 0.1 * meanTadSize) # using 1/10 * mean value for TAD size as cutoff for matching up and down peaks
      
      if(length(index)>0){
        upLoc = index%%length(upPeak[,1])
        downLoc = ceiling(index/length(upPeak[,1]))
        strap = data.frame(upPeak[upLoc, 1], downPeak[downLoc, 1])
        
        colMin = data.frame(apply(strap, 1, FUN=min)-5)
        colMin[,2] = 1
        colMax = data.frame(apply(strap, 1, FUN=max)+5)
        colMax[,2] = ncol(map)
        
        strap[,3] = apply(colMin, 1, FUN = max)
        strap[,4] = apply(colMax, 1, FUN = min)
        
        #strap = addFC(strap)
        colnames(strap) = c('up_loc', 'down_loc', 'leftEdge', 'rightEdge')
        strap[nrow(strap)+1,] <- NA
        return(list(strap, 0.1*meanTadSize))
      }
      else{
        strap = data.frame(matrix(ncol = 4, nrow = 0))
        colnames(strap) = c('up_loc', 'down_loc', 'leftEdge', 'rightEdge')
        strap[1,] <- NA
        return(list(strap, 1))
      }
      
    }
    else{
      strap = data.frame(matrix(ncol = 4, nrow = 0))
      colnames(strap) = c('up_loc', 'down_loc', 'leftEdge', 'rightEdge')
      strap[1,] <- NA
      return(list(strap, 1))
    }
  }
  
  reduceStrap <- function(strap, matchCutoff){
    if(nrow(strap)>2){
      strapReduce = data.frame(matrix(ncol = 4, nrow = 0))
      colnames(strapReduce) = c('up_loc', 'down_loc', 'leftEdge', 'rightEdge')
      
      flag = FALSE
      strapTmp = data.frame(matrix(ncol = 4, nrow = 0))
      colnames(strapTmp) = c('up_loc', 'down_loc', 'leftEdge', 'rightEdge')
      strapTmp[1, ] = strap[1, ]
      
      for(i in 2:(nrow(strap)-1)){
        if(abs(strap[i,1]-strap[(i-1),2]) < matchCutoff){
          if(flag){
            strapTmp[1,] <- c(strapTmp[1,1], strap[i,2], strapTmp[1,3], strap[i,4])
            colnames(strapTmp) = c('up_loc', 'down_loc', 'leftEdge', 'rightEdge')
          }
          else{
            strapTmp[1,] <- c(strap[(i-1),1], strap[i,2], strap[(i-1),3], strap[i,4])
            colnames(strapTmp) = c('up_loc', 'down_loc', 'leftEdge', 'rightEdge')
          }
          flag = TRUE
        }
        else{
          strapReduce = rbind(strapReduce, strapTmp)
          strapTmp = data.frame(matrix(ncol = 4, nrow = 0))
          colnames(strapTmp) = c('up_loc', 'down_loc', 'leftEdge', 'rightEdge')
          strapTmp[1, ] = strap[i, ]
          flag = FALSE
        }
      }
      strapReduce = rbind(strapReduce, strapTmp)
      strapReduce[nrow(strapReduce)+1,]<-NA
      return(strapReduce)
    }
    else{
      return(strap)
    }
  }
  
  bcpLeft = strapCalling(contactMap, upLimit, outputPath, outputName1, outputName2, startPoint, endPoint, window)
  oneSideLeft = oneSideData(bcpLeft)
  upSideLeft = oneSideLeft[[1]]
  downSideLeft = oneSideLeft[[2]]
  matchPeaks.temp = matchPeaks(upSideLeft, downSideLeft, p, contactMap)
  strapLeft = matchPeaks.temp[[1]]
  matchCutoff = matchPeaks.temp[[2]]
  
  strapReduce = reduceStrap(strapLeft, matchCutoff)
  
  if(outputName3!='NA'){
    pdf(outputName3,width = 4.5, height = 3.5)
    plot(upSideLeft)
    dev.off()
  }
  
  if(outputName4!='NA'){
    pdf(outputName4,width = 4.5, height = 3.5)
    plot(downSideLeft)
    dev.off()
  }
  return(list(strapReduce, upSideLeft, downSideLeft))
}

compareStrap <- function(list1, list2){
  # list1: a list for the first strap generating by getStrap function
  # list2: a list for the second strap generating by getStrap function
  # output: a list of result for the first and second strap with differential p value
  
  foldChange <- function(strap, upData, downData){
    # input: strap: strap of the first sample
    # upData: upData for the second sample
    # downData: downData for the second sample
    strap2 = head(strap, -1)
    getUp <- function(location){
      return(max(upData$posterior.mean[location[3]:location[4]]))
    }
    getDown <- function(location){
      return(min(downData$posterior.mean[location[3]:location[4]]))
    }
    
    getPvalue <- function(location){
      return(max(upData$posterior.prob[location[3]:location[4]]))
    }
    
    strap2[,ncol(strap2) + 1] = apply(strap2, 1, FUN=getUp)
    strap2[,ncol(strap2) + 1] =  apply(strap2, 1, FUN=getDown)
    strap2[,ncol(strap2) + 1] = log(abs(strap2[,ncol(strap2)-1])/abs(strap2[,ncol(strap2)]))
    strap2[,ncol(strap2) + 1] =  apply(strap2, 1, FUN=getPvalue)
    strap2[(nrow(strap2)+1),] <- NA
    return((strap2))
  }
  
  strap1 = foldChange(list1[[1]], list1[[2]], list1[[3]])
  strap1 = foldChange(strap1, list2[[2]], list2[[3]])
  strap2 = foldChange(list2[[1]], list2[[2]], list2[[3]])
  strap2 = foldChange(strap2, list1[[2]], list1[[3]])
  
  strap1$pvalue = pnorm(strap1[,7]-strap1[,11],0,2*0.5538227^2)
  colnames(strap1) <- c('upPeak.loc', 'downPeak.loc', 'leftEdge', 'rightEdge', 'upPeak.sample1', 
                        'downPeak.sample1', 'logFoldChange.sample1', 'strap.pValue.sample1', 'upPeak.sample2', 
                        'downPeak.sample2', 'logFoldChange.sample2', 'strap.pValue.sample2', 'diffStrap.pValue')
  
  strap2$pvalue = pnorm(strap2[,7]-strap2[,11],0,2*0.5538227^2)
  colnames(strap2) <- c('upPeak.loc', 'downPeak.loc', 'leftEdge', 'rightEdge', 'upPeak.sample2', 
                        'downPeak.sample2', 'logFoldChange.sample2', 'strap.pValue.sample2', 'upPeak.sample1', 
                        'downPeak.sample1', 'logFoldChange.sample1', 'strap.pValue.sample1', 'diffStrap.pValue')
  return(list(strap1, strap2))
}

diffStrap <- function(contactMap1, contactMap2, upLimit, window){
  # contactMap1: the first matrix of the Hi-C contact Map
  # contactMap2: the second matrix of the Hi-C contact Map
  # upLimit: up limitation of contact map
  # output: 
  strap1 = getStrap(contactMap1, upLimit, window)
  strap2 = getStrap(contactMap2, upLimit, window)
  compare1 = compareStrap(strap1, strap2)
  compare1[[1]]$direction = 'left'
  compare1[[2]]$direction = 'left'
  
  contactMap1.reverse = t(contactMap1[nrow(contactMap1):1, ncol(contactMap1):1])
  contactMap2.reverse = t(contactMap2[nrow(contactMap2):1, ncol(contactMap2):1])
  strap1 = getStrap(contactMap=contactMap1.reverse, upLimit=20, window)
  strap2 = getStrap(contactMap=contactMap2.reverse, upLimit=20, window)
  compare2 = compareStrap(strap1, strap2)
  compare2[[1]][,1:4] = nrow(contactMap1) - compare2[[1]][,1:4]
  compare2[[2]][,1:4] = nrow(contactMap2) - compare2[[2]][,1:4]
  compare2[[1]]$direction = 'right'
  compare2[[2]]$direction = 'right'
  
  pair1ToPair2 = rbind(compare1[[1]], compare2[[1]])
  pair1ToPair2 = pair1ToPair2[order(pair1ToPair2$upPeak.loc),]
  
  pair2ToPair1 = rbind(compare1[[2]], compare2[[2]])
  pair2ToPair1 = pair2ToPair1[order(pair2ToPair1$upPeak.loc),]
  
  return(list(head(pair1ToPair2, -1), head(pair2ToPair1, -1)))
}
