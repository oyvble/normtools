#test_preprocessENmix
#TESTING DIFFERENT IMPLEMENTATIONS OF EM_ESIMATE:


test_preprocessENmix = function() {
  require(ENmix)
  require(minfiData)
  require(normtools)
  
  #PRE-EXAMPLE:
  #Extract a sample
  set.seed(1)
  n=62000
  x = rgamma(n,100,1/1000)
  #COMPARE TIME ESIMATES
  library(microbenchmark) #compare speed
  microbenchmark(
    ENmix:::EM_estimate(x), #390
    normtools::EM_estimateR(x), #172 (2.3x)
    normtools::EM_estimateC(x), #52 (7.5x)
    times=10
  )
  
  #TWST METHOLDS  
  nCores = detectCores() #number of cores to utilize
  #nCores=4
  #READ DATA:
  path <- file.path(find.package("minfiData"),"extdata")
  rgSet <- readidat(path = path,recursive = TRUE)
  #rgSet <- rgSet[,1:4]
  
  #TEST RE-IMPLEMENTATION OF preprocessENmix
  bgMethods = c("oob","est","neg") #traverse each method
  timeList <- valList <- list()
  for(bgParaEst in bgMethods) {
#bgParaEst="oob"
    time = system.time({
      mdat <-ENmix::preprocessENmix(rgSet, bgParaEst=bgParaEst, dyeCorr="none", nCores=nCores)
    })[3]
    timeR = system.time({
      mdatR <-normtools::preprocessENmixR(rgSet, bgParaEst=bgParaEst, dyeCorr="none", nCores=nCores)
    })[3]
    timeC = system.time({
      mdatC <-normtools::preprocessENmixC(rgSet, bgParaEst=bgParaEst, dyeCorr="none", nCores=nCores)
    })[3]
    timeList[[bgParaEst]] = c(time,timeR,timeC) #(155, 76, 55)
    
    #COMPARE:
    mean0 = mean(assays(mdat)$Meth+assays(mdat)$Unmeth)
    meanR = mean(assays(mdatR)$Meth+assays(mdatR)$Unmeth)
    meanC = mean(assays(mdatC)$Meth+assays(mdatC)$Unmeth)
    
    valList[[bgParaEst]] = c(mean0,meanR,meanC)
  }
}