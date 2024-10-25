#THIS SCRIPT WILL FOR EACH IDAT FILE IN THE IDAT FOLDER DO FOLLOWING: 
#1) Copy file, 2) Read idat info (using illuminaio::readIDAT) 3) Rewrite copied file with replaced MEAN/SD as 0 for Specific probe type (defualt is SNP probes)
#AUTHOR: Oyvind Bleka 
#rm(list=ls());gc()
#setwd("C:\\Users\\oyvbl\\Documents\\Methylation\\rmSNPinIDAT") #SELECT WORKING DIRECTORY
#source("IDATreplacer.R")

#IDATfold = "IDATgz" #Folder name with IDAT files (possibly gz compressed)
#IDATfold2 = "IDAT_OUT" #Folder name with modidied IDAT files
#prefix = "SNPrm_" #prefix name (to name copied files). CAN ALSO BE "" (empty)
#compressToGz = TRUE #TRUE/FALSE SHOULD THE MODIFIED FILES BE COMPRESSED TO Gz again?

IDATreplacer = function(IDATfold,IDATfold2, probeType = "SNP", prefix= paste0(probeType,"rm_"),compressToGz = TRUE, verbose=TRUE) {
  library(illuminaio) #USE illuminaio R package to read idat files
  #if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  #if(!require(illuminaio)) BiocManager::install("illuminaio", version = "3.8")
  if(!file.exists(IDATfold2))  dir.create(IDATfold2)
  
  #OBTAIN ADDRESS OF WHAT TO REMOVE (MEAN,SD set to zero)
  snpaddfile =  paste0( path.package("normtools"),"/SNPaddressinfo.csv") # Get file name in package 
  #snpaddfile="SNPaddressinfo.csv" #snp address file (lookup table)
  #following what = "all"
  
  #GET ADRESS OF SNPS:
  addr = read.table(snpaddfile,head=TRUE,sep=",",as.is=TRUE)
  addr = as.numeric(unlist(strsplit(addr[,3],"/"))) #get adresses
  
  
  rgSet = minfi::read.metharray(ff,verbose=TRUE,extended=TRUE)
  dat <- minfi::getProbeInfo(rgSet, type = typ )
  
  #Get info about snps
  tab = numeric()
  SNPtypes = c("SnpI","SnpII")
  for(typ in SNPtypes) {
  #  typ= SNPtypes[1]
    dat <- getProbeInfo(rgSet, type = typ)
    if( typ == "SnpI") add <- cbind(dat[,1], paste0( dat[,2],"/",dat[,3] ))  #get addresses
    if( typ == "SnpII") add <- dat[,1:2]  #get addresses
    tab = rbind(tab, cbind(typ,as.matrix(add)))
  }
  colnames(tab) = c("SNPtype","Site","Address")
  
  
  #GET idat FILES
  ff <- list.files(IDATfold) #get file names
  
  #FOR EACH FILE IN IDATFOLD
  for(filename in ff) { 
    # file=ff[1]
    if(verbose) print(paste0("PROCESSING FILE ",filename))
    ofile <- paste0(prefix,filename) #this is copied file (to SNP values to zero)
    file2 = paste0(IDATfold,"\\",filename) #with path
    ofile2 = paste0(IDATfold2,"\\",ofile)#with path
    file.copy(file2,ofile2) #COPY FILE FIRST (MAKE A DUPLICATE):
    obj = illuminaio::readIDAT(file2)
  
    #all(addr%in%obj$MidBlock) #all found ins elem? 
    #all(rownames(obj$Quants)==obj$MidBlock) #same order?
    Mean = obj$Quants[,1] #get Mean reads
    SD = obj$Quants[,2] #get SD reads
  
    #RESET SNP INFO:
    rmrow =obj$MidBlock%in%addr #get index of SNPs
    Mean[rmrow] = 0 #SET MEAN SIGNAL TO ZERO
    SD[rmrow] = 0 #SET SD SIGNAL TO ZERO (also important)
  
    #REPLACE INFO in copied file: 
    if(grepl(".gz",ofile2)) { #check if file must be unzipped
    	require(R.utils) #install.packages("R.utils")
    	gunzip(ofile2)
    	ofile2 = gsub(".gz","",ofile2) #update file name
    }
    con <- file(ofile2, "r+b") #both read and write in binary mode
  
    #GET POSITION TO CHANGE
    where <- obj$fields["Mean", "byteOffset"]
    seek(con, where = where, origin = "start",rw="write")
    writeBin(as.integer(Mean),con, size=2, endian="little")  #Mean vector must be integer
  
    where <- obj$fields["SD", "byteOffset"]
    seek(con, where = where, origin = "start",rw="write")
    writeBin(as.integer(SD),con, size=2, endian="little")  #SD vector must be integer
    close(con)
  
    if(compressToGz) {
    	require(R.utils) #install.packages("R.utils")
    	gzip(ofile2) #Compress to gz
    }
  }  #End for each file
} 

