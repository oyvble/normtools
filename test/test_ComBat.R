
library(sva)
library(bladderbatch)
data(bladderdata)
dat <- bladderEset[1:50,]

pheno = pData(dat)
edata = exprs(dat)
batch = pheno$batch
mod = model.matrix(~as.factor(cancer), data=pheno)
ref.batch = 1

#dat=edata;mod=NULL;par.prior=TRUE;mean.only = FALSE;maxGB=2; log2transform=FALSE; transposed=FALSE

#NO REFERENCE BATCH

#parametric adjustment (no covariates)
par1a = sva::ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE)
par1b = normtools::ComBat2(dat=edata, batch=batch, mod=NULL, par.prior=TRUE)
print(min(abs(par1a-par1b)))

#parametric adjustment (include factor in covariates)
par2a = sva::ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE)
par2b = normtools::ComBat2(dat=edata, batch=batch, mod=mod, par.prior=TRUE)
print(min(abs(par2a-par2b)))

# non-parametric adjustment, (no covariates)
np1a = sva::ComBat(dat=edata, batch=batch, mod=NULL, par.prior=FALSE)
np1b = normtools::ComBat2(dat=edata, batch=batch, mod=NULL, par.prior=FALSE)
print(min(abs(np1a-np1b)))

#parametric adjustment (include factor in covariates)
np2a = sva::ComBat(dat=edata, batch=batch, mod=mod, par.prior=FALSE)
np2b = normtools::ComBat2(dat=edata, batch=batch, mod=mod, par.prior=FALSE)
print(min(abs(np2a-np2b)))

#SET BATCH 1 AS REFERENCE BATCH

#parametric adjustment (no covariates)
par1a = sva::ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE,ref.batch = ref.batch)
par1b = normtools::ComBat2(dat=edata, batch=batch, mod=NULL, par.prior=TRUE,ref.batch = ref.batch)
print(min(abs(par1a-par1b)))

#parametric adjustment (include factor in covariates)
par2a = sva::ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE,ref.batch = ref.batch)
par2b = normtools::ComBat2(dat=edata, batch=batch, mod=mod, par.prior=TRUE,ref.batch = ref.batch)
print(min(abs(par2a-par2b)))

# non-parametric adjustment, (no covariates)
np1a = sva::ComBat(dat=edata, batch=batch, mod=NULL, par.prior=FALSE,ref.batch = ref.batch)
np1b = normtools::ComBat2(dat=edata, batch=batch, mod=NULL, par.prior=FALSE,ref.batch = ref.batch)
print(min(abs(np1a-np1b)))

#parametric adjustment (include factor in covariates)
np2a = sva::ComBat(dat=edata, batch=batch, mod=mod, par.prior=FALSE,ref.batch = ref.batch)
np2b = normtools::ComBat2(dat=edata, batch=batch, mod=mod, par.prior=FALSE,ref.batch = ref.batch)
print(min(abs(np2a-np2b)))



