
#Estimating number of modes per marker by traversing a grid (vectorized)
#Idea: Calculate percentage of points in each segment of grid.
#Iterate each segment sequentically and store maximum percentages until a segment falls within a percentage (then a mode has been identified)
estNmodes = function(X,req = 0.1,bw=0.1, minCount=1,transposed=FALSE, verbose=TRUE) {
  #req Criterion of %-drop from maximum-proportion to count as mode
  #bw band-width of grid
  #minCount minimum number of samples on a bin to accept
  #transposed whether the input matrix is transposed or not
 
  breaks = seq(0,1,bw)
  nBreaks = length(breaks) #number of breaks
  
  if(transposed) {  #if X is (samples x sites)
    nData = nrow(X) #number of data per marker
    nSites = ncol(X) #number of sites
  } else {
    nData = ncol(X) #number of data per marker
    nSites = nrow(X) #number of sites
  }
  #Iterate each segment
  nModes = rep(0,nSites) #number of sites
  hasCounted = rep(FALSE,nSites) #whether a mode has been counted after obtaining new maximum
  isBelow =  rep(FALSE,nSites)  #default info about whether is below threshold
  
  for(i in 2:length(breaks)) {
    x0 = breaks[i-1]
    x1 = breaks[i]
    bool = X>=x0 & X<x1 #obtain true/false matrix
    if(transposed) {
      counts = colSums(bool) #insert
    } else {
      counts = rowSums(bool) #insert
    }
    
    #preproposess:
    counts[counts<=minCount] = 0 #set count to zero when below minCount
    
    prop = counts/nData #proportion of data within segment
    
    if(i>2) { #UPDATING MAX IF CLIMBING
      #Identify Hill-climbing situations:
      isGreater = prop > maxVal #is greater than 
      maxVal[isGreater] = prop[isGreater] #insert
      
      #Reset counter when doing new climb
      hasCounted[isGreater] = FALSE
      
      #Check if previous was below thresh
      prevBelow = isBelow #check if prev was also below threshold
      
      #Identify Hill-descendant situations:
      isBelow = prop < req*maxVal #whether dropping below threshold
      
      #identify which sites should update mode count: 
      #Current is below AND Previous was not below AND it has not been counted earlier
      updateInd = !prevBelow & isBelow & !hasCounted 
      
      nModes[updateInd] = nModes[updateInd] + 1 #count mode when dropping below thershold
      hasCounted[updateInd] = TRUE #counted
      maxVal[updateInd] = prop[updateInd] #set new max if below
      
    } else {  #init maxval to current prop
      maxVal = prop #store max value
    }
    if(verbose) print(paste0( round(i/length(breaks)*100),"% complete"))
    
  }
  nModes[!hasCounted] = nModes[!hasCounted] + 1 #count mode if not earlier counted (because of reaching boarder)
  return(nModes)  
}

