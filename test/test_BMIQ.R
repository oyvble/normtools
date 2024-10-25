#SCRIPT FOR COMPARING BMIQ implementations
#rm(list=ls())
require(minfiData)
library(normtools)

#Must define seed since bmiq is based on drawing random markers
seed=1

#READ DATA:
fn <- list.files(file.path(find.package("minfiData"),"extdata"),recursive = T,full=T)
rgSet = minfi::read.metharray( fn[grepl("Grn.idat",fn)]  ,verbose=TRUE,extended=TRUE) #Read IDAT files to R object. "Extended" allows quality control
#rgSet <- readidat(path = path,recursive = TRUE)
mdat <- preprocessNoob(rgSet) #obtaining methylated data (avoids nan)
bdat <- getBeta(mdat) #obtaining methylated data
siteRM = unique(which(rowSums(is.na(bdat))>0)) #remove following sites because of nan 
if(length(siteRM)) bdat <- bdat[siteRM,] #remove bad sites

#Performing BMIQ
set.seed(seed)
time1 = system.time({ bdatBMIQ1 = normtools::bmiq.mc(bdat,implementation = "watermelon")})[3] #120s
set.seed(seed)
time2 = system.time({ bdatBMIQ2 = normtools::bmiq.mc(bdat,implementation = "fast")})[3] #7s (17x faster)

#OUTPUT COMPARISON
siteprint = 135476:135500 #first is Type1, remaining is Type2
i=1
cbind(bdat[siteprint,i],bdatBMIQ1[siteprint,i],bdatBMIQ2[siteprint,i])
diff = bdatBMIQ1[siteprint,i]-bdatBMIQ2[siteprint,i]

#class(rgSet) is rgDataSet (go to corresponding block in preprocessENmix)
#COMPARING PART OF CODE IN preprocessENmix:

library(IlluminaHumanMethylation450kmanifest)
manifestType = setNames(Manifest$Type,Manifest$Name) #obtain probe type
indType2 = rownames(bdat)%in%names(manifestType[manifestType=="II"])
plot(density(bdat[indType2,i]),main=colnames(bdat)[i])
lines(density(bdatBMIQ1[indType2,i]),col=2)
lines(density(bdatBMIQ2[indType2,i]),col=3)


