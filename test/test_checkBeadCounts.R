library(ewastools)
library(normtools)

idat_files = paste0("C:\\Users\\oyvbl\\Documents\\Methylation\\GSE87650\\data\\IDAT\\",
c("GSM2336819_9647455142_R02C01","GSM2336827_9647455142_R04C02","GSM2336833_9647450103_R04C01","GSM2336856_9647455070_R03C01"))

meth = ewastools::read_idats(idat_files,quiet=TRUE) # `quiet=TRUE` supresses the progress bar
failed = checkBeadCounts(meth)
