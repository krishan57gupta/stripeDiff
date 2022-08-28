################################################################ simulation code
############### define simulation function ##################
simulateStrap <- function(l, m, slope, strap.hight, margin.strap = 0.10, margin.boudary= 0.10, 
                          strapType = 'left', matrixType = 'noLoop', strapStart = 1, strapEnd = l){
  # l: length of TAD
  # m: mean of interaction T1: 3 ~ 4 T2: 1~2
  # slope: log(interaction) VS log(distance) -0.65 ~ -0.85
  # strap.hight: hight of the strap (>=1), set to 1 as not strap
  # strapType: set "left", "right", "None" as left strap, right strap, or none strap (need to add loop type)
  # matrixType: "loop" or "noLoop"
  # return: simulation TAD, mask for strap, mask for boudary
  require(sigmoid)
  require(geosphere)
  
  sgmd <- function(x, n, direction){
    if(direction == "strap"){
      y = logistic(x/n*10, x0 = 5, k=1)
    }
    else{
      y = logistic(x/n*10, x0 = 5, k=1)
    }
    return(y)
  }
  
  strapMask <- function(strapType, l, strap1, strap2, strap.hight, strapStart, strapEnd){
    mask1 = matrix(1, nrow=l, ncol=l)
    mask2 = matrix(1, nrow=l, ncol=l)
    
    if(strapType == 'None'){
      # mask for strap
      for(i in strapStart:l){
        for(j in 1:strap1){
          if(i!=j){
            h2 = sgmd(j, strap1, 'strap')
            mask1[i, j] =  h2
          }
          else{
            h2 = sgmd(j, strap1, 'strap')
            mask1[i, j] = h2
          }
        }
      }
      
      # mask for boundary
      for(i in 1:l){
        for(j in l:(l-strap2+1)){
          if(i!=j){
            h2 = sgmd((l-j+1), strap2, 'boundary')
            mask2[i, j] =  h2
          }
          else{
            h2 = sgmd((l-j+1), strap2, 'boundary')
            mask2[i, j] = h2
          }
        }
      }
      
    }
    else{
      if(strapType == 'left'){
        # mask of left strap for TAD
        # mask for strap
        for(i in strapStart:l){
          for(j in 1:strap1){
            if(i!=j){
              h2 = sgmd(j, strap1, 'strap')
              mask1[i, j] =  h2*strap.hight
            }
            else{
              h2 = sgmd(j, strap1, 'strap')
              mask1[i, j] = h2*strap.hight
            }
          }
        }
        
        # mask for boundary
        for(i in 1:l){
          for(j in l:(l-strap2+1)){
            if(i!=j){
              h2 = sgmd((l-j+1), strap2, 'boundary')
              mask2[i, j] =  h2
            }
            else{
              h2 = sgmd((l-j+1), strap2, 'boundary')
              mask2[i, j] = h2
            }
          }
        }
      } 
      else{
        # mask of right strap for TAD
        for(i in 1:strapEnd){
          for(j in 1:strap1){
            if(i!=j){
              h2 = sgmd(j, strap1, 'strap')
              mask1[i, j] =  h2*strap.hight
            }
            else{
              h2 = sgmd(j, strap1, 'strap')
              mask1[i, j] = h2*strap.hight
            }
          }
        }
        
        # mask for boundary
        for(i in 1:l){
          for(j in l:(l-strap2+1)){
            if(i!=j){
              h2 = sgmd((l-j+1), strap2, 'boundary')
              mask2[i, j] =  h2
            }
            else{
              h2 = sgmd((l-j+1), strap2, 'boundary')
              mask2[i, j] = h2
            }
          }
        }
        
        mask1 = t(mask1[,ncol(mask1):1])
        mask2 = t(mask2[,ncol(mask2):1])
      }
    }
    return(list(mask1, mask2))
  }
  
  simulateTAD <- function(l, type, strap){
    # simulate two type of TADs, with loop or without loop
    # l: size of TAD
    # type: type of TADs, "loop" or "noLoop"
    # return: TAD matrix
    
    getTad <- function(l){
      TAD =  matrix(0, nrow = l, ncol = l)
      for(i in 1:l){
        for(j in 1:i){
          if(i!=j){
            m1 = log(abs(i - j))*(slope) + m
            m2 = log(abs(i - j))*(0.17) + 0.15
            r = rnorm(1, mean = m1, sd = m2^2)  # original sd = m2^2
            TAD[i, j] = TAD[i, j] + exp(1)^(r)
          }
          else{
            m1 = m
            m2 = 0.15
            r = rnorm(1, mean = m1, sd = m2^2)  # original sd = m2^2
            TAD[i, j] = TAD[i, j] + exp(1)^(r)
          }
        }
      }
      return(TAD)
    }
    
    getLoop <- function(l, strap){
      TAD =  matrix(0, nrow = l, ncol = l)
      for(i in 1:l){
        for(j in 1:i){
          m1 = log(abs(l-i + j))*(slope/2) + m
          m2 = log(abs(l-i + j))*(0.17) + 0.15
          r = rnorm(1, mean = m1, sd = m2^2)  # original sd = m2^2
          
          if( ((i-l + strap)^2 + (j - strap)^2)^0.5 <= strap ){
            TAD[i, j] = TAD[i, j] + exp(1)^(r)*5 
          }
          else{
            TAD[i, j] = TAD[i, j] + exp(1)^(r)
          }
        }
      }
      
      TAD[upper.tri(TAD)] = t(TAD)[upper.tri(t(TAD))]
      return(TAD)
    }
    
    if(type == 'noLoop'){
      TAD0 = getTad(l)
      TAD = TAD0
      loopTAD = TAD
    }
    else{
      ##
      TAD0 = getTad(l)
      loopTAD = getLoop(l, strap)
      TAD = (1-0.5)*TAD0 + 0.5*loopTAD
    }
    
    return(list(TAD, loopTAD, TAD0))
  }
  
  strap1 = ceiling(l*margin.strap)
  strap2 = ceiling(l*margin.boudary)
  TAD.temp = simulateTAD(l, matrixType, 2)
  TAD = TAD.temp[[1]]
  
  mask1 = strapMask(strapType, l, strap1, strap2, strap.hight, strapStart, strapEnd)[[1]]
  mask2 = strapMask(strapType, l, strap1, strap2, strap.hight, strapStart, strapEnd)[[2]]
  
  TAD = TAD*mask1*t(mask2)
  TAD[upper.tri(TAD)] = t(TAD)[upper.tri(t(TAD))]
  return(list(TAD, mask1, mask2, TAD.temp[[2]], TAD.temp[[3]]))
}

simulateBackground <- function(l, m, slope, sd, segment, jointType, tier2Hight, strapType, strapHight, 
                               strapWidth, loopType){
  # l: length of TAD
  # m: mean of interaction T1: 3 ~ 4 T2: 1~2
  # slope: log(interaction) VS log(distance) -0.65 ~ -0.85
  # segment: a vector of length for each TAD in the background
  # jointType: a vector for whether jointing two TADs nearby each other
  # tier2Hight: range of hight for tier 2 TAD (a vecter)
  # strapType: a vector for the strap type of each TAD ('left', 'right', 'None')
  # strapHight: a vector for the strap hight ( 2 ~ 2.5)
  # strapWidth: a vector for the strap width (0.03 ~ 0.06)
  # loopType: a vector for the loop type ('noLoop', 'loop')
  # return: simulated background TAD
  
  TAD =  matrix(0, nrow = l, ncol = l)
  for(i in 1:l){
    for(j in 1:l){
      if(i!=j){
        m1 = log(abs(i - j))*(slope) + m
        m2 = log(abs(i - j))*(0.17) + sd               
        r = rnorm(1, mean = m1, sd = m2^2) # original sd = m2^2
        TAD[i, j] = TAD[i, j] + exp(1)^(r)
      }
      else{
        m1 = m
        m2 = sd               
        r = rnorm(1, mean = m1, sd = m2^2)  # original sd = m2^2
        TAD[i, j] = TAD[i, j] + exp(1)^(r)
      }
    }
  }
  
  for(i in 1:(length(segment)-1)){
    if(jointType[i] == 'yes'){
      size = segment[i] + segment[i+1]
      if(strapType[i] == 'left'){
        # for the left strap
        TAD.temp = simulateStrap(size, 1.5, -0.7, strapHight[i], margin.strap = strapWidth[i], margin.boudary= 0.05, 
                                 strapType = strapType[i], matrixType = loopType[i], strapStart = segment[i])[[1]]
      }
      else{
        # for the right strap
        TAD.temp = simulateStrap(size, 1.5, -0.7, strapHight[i], margin.strap = strapWidth[i], margin.boudary= 0.05, 
                                 strapType = strapType[i], matrixType = loopType[i], strapEnd = segment[i])[[1]]
      }
      h = runif(1, tier2Hight[1], tier2Hight[2])
      TAD[(sum(segment[0:(i-1)])+1):sum(segment[0:(i+1)]), (sum(segment[0:(i-1)])+1):sum(segment[0:(i+1)])] = 
        TAD.temp*h + TAD[(sum(segment[0:(i-1)])+1):sum(segment[0:(i+1)]), (sum(segment[0:(i-1)])+1):sum(segment[0:(i+1)])]
    }
  }
  
  TAD[upper.tri(TAD)] = t(TAD)[upper.tri(t(TAD))]
  return(TAD)
}

simulateContactMatrix <- function(backGroudMap, l, h, S1, S2, m, d, strapType, tadType){
  # backGroudMap: a contact matrix of back ground
  # l: a vector of TAD size
  # h: a vector of fold change between strap contact and TAD contact
  # S1: a vector of strap width
  # S2: a vector of boundary width
  # m: a vector of mean of interaction T1: 3 ~ 4 T2: 1~2
  # d: a vector of log(interaction) VS log(distance) -0.65 ~ -0.85
  # strapType: a vector for type of each strap
  # tadType: a vector for type of each TAD
  # return: simulated contact matrix
  
  TAD.temp = backGroudMap
  TAD0 = TAD.temp
  start = 1
  end = l[1]
  
  for(i in 1:length(l)){
    TAD = simulateStrap(l[i], m[i], d[i], h[i], margin.strap = S1[i], margin.boudary = S2[i], strapType = strapType[i], matrixType = tadType[i])[[1]]
    # m1 = mean(TAD0[start:end, 1])
    # TAD0[start:end, start:end] = TAD + m1
    TAD0[start:end, start:end] = TAD + TAD0[start:end, start:end]
    if(i!=length(l)){
      start = end+1
      end = start+l[i+1]-1
    }
  }
  
  return(TAD0)
}

generateParameter <- function(){
  l = c()
  h = c()
  S1 = c()
  S2 = c()
  m = c()
  d = c()
  strapType = c()
  matrixType = c()
  
  for(i in 1:ntad){
    l = c(l, ceiling(runif(1, ft, lt)))
    h = c(h, runif(1, 1.8, 2.2))
    S1 = c(S1, runif(1, 0.05, 0.12))
    S2 = c(S2, runif(1, 0.05, 0.1))
    m = c(m, runif(1, 1, 3))
    d = c(d, -runif(1, 0.65, 0.85))
    
    seed = runif(1, 0, 1)
    if(seed< 0.3){
      strapType = c(strapType, 'None')
      if(seed< 0.15){
        matrixType = c(matrixType, 'loop')
      }else{
        matrixType = c(matrixType, 'noLoop')
      }
    }else if(seed < 0.65){
      strapType = c(strapType, 'left')
      if(seed< 0.47){
        matrixType = c(matrixType, 'loop')
      }else{
        matrixType = c(matrixType, 'noLoop')
      }
    }else{
      strapType = c(strapType, 'right')
      if(seed< 0.82){
        matrixType = c(matrixType, 'loop')
      }else{
        matrixType = c(matrixType, 'noLoop')
      }    
    }
  }
  
  jointType = c()  # yes, no
  strapType.backGround = c() # left, right, None
  loopType.backGround = c() # loop, noLoop
  strapHight.backGround = c() # 1, 2 ~ 2.5
  strapWidth.backGround = c()
  
  for(i in 1:ntad){
    seed = runif(1, 0, 1)
    if(seed < 0.3){
      jointType = c(jointType, 'yes')
      
      if(i != ntad){
        if(strapType[i] == 'None'){
          if(runif(1, 0, 1) < 0.4){
            strapType.backGround = c(strapType.backGround, 'left')
          } else if(runif(1, 0, 1) < 0.4 & runif(1, 0, 1) > 0.2){
            strapType.backGround = c(strapType.backGround, 'right')
          }else{
            strapType.backGround = c(strapType.backGround, 'None')
          }
        }else if(strapType[i] == 'left'){
          if(runif(1, 0, 1) < 0.3){
            strapType.backGround = c(strapType.backGround, 'left')
          }else{
            strapType.backGround = c(strapType.backGround, 'None')
          }
        }else{
          if(runif(1, 0, 1) < 0.3){
            strapType.backGround = c(strapType.backGround, 'right')
          }else{
            strapType.backGround = c(strapType.backGround, 'None')
          }
        }
      }
      else{
        strapType.backGround = c(strapType.backGround, 'None')
      }
      
      loopType.backGround = c(loopType.backGround, 'loop')
      strapHight.backGround = c(strapHight.backGround, runif(1, 2, 2.5))
      strapWidth.backGround = c(strapWidth.backGround, runif(1, 0.03, 0.06))
    }else{
      jointType = c(jointType, 'no')
      strapType.backGround = c(strapType.backGround, 'None')
      loopType.backGround = c(loopType.backGround, 'noLoop')
      strapHight.backGround = c(strapHight.backGround, 1)
      strapWidth.backGround = c(strapWidth.backGround, 0.06)
    }
  }
  
  parameter = data.frame(l, h, S1, S2, m, d, strapType, matrixType, jointType, strapType.backGround, 
                         loopType.backGround, strapHight.backGround, strapWidth.backGround)
  return(parameter)
}

pairParameter <- function(parameter){
  changePoint = data.frame(bin, 'None')
  colnames(changePoint) <- c('loc', 'type')
  nonChange = data.frame(bin, 'None')
  colnames(nonChange) <- c('loc', 'type')
  position = c(1, cumsum(parameter$l))
  for(i in 1:nrow(parameter)){
    seedLeft = runif(1, 0, 1)
    seedRight = runif(1, 0, 1)
    if(parameter$strapType[i] == 'left' ){
      if(seedLeft < 0.5){
        parameter$strapType[i] = 'None'
        
        changePoint.temp = data.frame(position[i], 'left')
        colnames(changePoint.temp) <- c('loc', 'type')
        
        changePoint = rbind(changePoint, changePoint.temp)
      }else{
        nonChange.temp = data.frame(position[i], 'left')
        colnames(nonChange.temp) <- c('loc', 'type')
        
        nonChange = rbind(nonChange, nonChange.temp)
      }
      
    }
    
    if(parameter$strapType[i] == 'right'){
      if(seedRight < 0.5){
        parameter$strapType[i] = 'None'
        
        changePoint.temp = data.frame(position[i+1], 'right')
        colnames(changePoint.temp) <- c('loc', 'type')
        
        changePoint = rbind(changePoint, changePoint.temp)
      }else{
        nonChange.temp = data.frame(position[i+1], 'right')
        colnames(nonChange.temp) <- c('loc', 'type')
        
        nonChange = rbind(nonChange, nonChange.temp)
      }
    }
  }
  return(list(parameter, changePoint[-1,], nonChange[-1,]))
}

############# generate simulation data ######################
generateSimulationData <-function(n1, n2, outputPath){
  for(i in n1:n2){
    print(i)
    ######  generate simulation data with strap (pair1) ###### 
    parameter = generateParameter()
    backGroudMap.0.0.1 = simulateBackground(sum(parameter$l), 1.5, -0.7, 0.15, parameter$l, 
                                            parameter$jointType, c(1, 1.2), parameter$strapType.backGround,
                                            parameter$strapHight.backGround, parameter$strapWidth.backGround, 
                                            parameter$loopType.backGround)
    map.0.0.1 = simulateContactMatrix(backGroudMap.0.0.1, parameter$l, parameter$h, parameter$S1, parameter$S2, 
                                      parameter$m, parameter$d, parameter$strapType, parameter$matrixType)
    ########################## increasing depth extra part done, originally not required
    aa=map.0.0.1
    N<-function(x,minv=0,maxv=1){
      (((x-min(x))/(max(x)-min(x)))*(maxv-minv))+minv
    }
    # aa=N(aa,0,sample(1000:100000,1))
    # aa=aa/sum(aa)*seqD*10^6
    map.0.0.1=aa
    ##########################
    write.table(map.0.0.1, file = paste(outputPath, '/simulation',i,'_pair1.txt',sep = ''),
                quote = F, row.names = F,col.names = F,sep = '\t')
    ###### generate simulation data without strap (pair2) ###### 
    changeStrap = pairParameter(parameter)
    parameter2 = changeStrap[[1]]
    backGroudMap.0.0.2 = simulateBackground(sum(parameter2$l), 1.5, -0.7, 0.15, parameter2$l, 
                                            parameter2$jointType, c(1, 1.2), parameter2$strapType.backGround,
                                            parameter2$strapHight.backGround, parameter2$strapWidth.backGround, 
                                            parameter2$loopType.backGround)
    map.0.0.2 = simulateContactMatrix(backGroudMap.0.0.2, parameter2$l, parameter2$h, parameter2$S1, parameter2$S2, 
                                      parameter2$m, parameter2$d, parameter2$strapType, parameter2$matrixType)
    ########################## increasing depth extra part done, originally not required
    aa=map.0.0.2
    N<-function(x,minv=0,maxv=1){
      (((x-min(x))/(max(x)-min(x)))*(maxv-minv))+minv
    }
    #aa=N(aa,0,sample(1000:100000,1))
    # aa=aa/sum(aa)*seqD*10^6
    map.0.0.2=aa
    ##########################
    write.table(map.0.0.2, file = paste(outputPath, '/simulation',i,'_pair2.txt',sep = ''),
                quote = F,row.names = F,col.names = F,sep = '\t')
    ###### write parameter and the location of changed strap ###### 
    write.table(parameter, file = paste(outputPath, '/simulation',i,'_parameter1.txt',sep = ''),
                quote = F, row.names = F,sep = '\t')
    write.table(parameter2, file = paste(outputPath, '/simulation',i,'_parameter2.txt',sep = ''),
                quote = F, row.names = F,sep = '\t')
    write.table(changeStrap[[2]], file =  paste(outputPath, '/simulation',i,'_changeStrap.txt',sep = ''),
                quote = F,row.names = F,col.names = T,sep = '\t')
    write.table(changeStrap[[3]], file =  paste(outputPath, '/simulation',i,'_nonChangeStrap.txt',sep = ''),
                quote = F,row.names = F,col.names = T,sep = '\t')
  }
}

suppressMessages(library("optparse"))
options(warn=-1)

option_list = list(
  make_option(c("-a", "--number1"), type="character",
              help="Start number of simulation files", metavar="character"),
  make_option(c("-b", "--number2"), type="character",
              help="End number of simulation files", metavar="character"),
  make_option(c("-o", "--outputFile"), type="character",
              help="Output path", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
generateSimulationData(opt$number1, opt$number2, opt$outputFile)

### running example
### cd /archive/tmhgxw22/Hi-C_strip/data/simulation/bin
### Rscript simulation_stripe.R -a 51 -b 100 -o /simulation/data