#!/usr/bin/env R

# Author: Sean Maden
# Get 1,000 probes with low variance across 7 noncancer tissues,
# then study genome mapping patterns to make Table S7.
#
# This script takes as input the list of probe data by tissue (read 
# from the ".*data.tables.filt.*" files available in supplement).
# A 2-step low variance analysis is performed across tissues, then 
# the 1,000 probes with lowest variances across tissues for both 
# methods are selected. Finally, probes with low mean DNAm filter 
# are retained as the final recurrent low-variance probes.
# 
# To see the construction of Table S7, refer to the script `tableS7.R`.
# To see how genome mapping summaries are plotting for Figure S8, see 
# the script `figS8.R`.

library(minfiData)

#----------
# load data
#----------
# get minfi-formatted probe annotations
cga <- getAnnotation(get(data("MsetEx")))

# load lfilt containing probe data tables by tissue
pkgname <- "recountmethylationManuscriptSupplement"
nct7.dir <- system.file("extdata", "nct7", package = pkgname)
lf.fn <- "lfilt_cgtables_nct.rda"
lfilt <- get(load(file.path(nct7.dir, lf.fn)))

#---------------------------
# quantile mvps (bin method)
#---------------------------
# params
qiv <- seq(0, 1, 0.01) ; qwhich <- c(11)
binv <- seq(0, 1, 0.01)[1:100] # binned bval mean
lmvp <- list()
# get quantile filtered probes
for(t in 1:length(lfilt)){
  ba <- lfilt[[t]]; lci <- list(c())
  # iterate on betaval bins
  for(b in binv){
    bf <- ba[ba$mean >= b & ba$mean < (b + 0.01), ] # filtered on mean bin
    qf <- quantile(bf$var, qiv)[qwhich]
    # iterate on ci values
    for(q in 1:length(qf)){lci[[q]] <- c(lci[[q]], bf$cgid[bf$var < qf[q]])}
    message(t, ":", b)
  }
  names(lci) <- paste0("ci:", names(qf)); lmvp[[names(lfilt)[t]]] <- lci
}
lmvp.bin <- lmvp
#----------------------------------
# quantile mvps (abs cutoff method)
#----------------------------------
lmvp <- list()
qiv <- seq(0, 1, 0.01); qwhich <- c(11)
binv <- seq(0, 1, 0.01)[1:100] # binned bval mean
# get quantile filtered probes
for(t in 1:length(lfilt)){
  ba <- lfilt[[t]]; lci <- list(c())
  qf <- quantile(ba$var, qiv)[qwhich]
  for(q in 1:length(qf)){lci[[q]] <- ba$cgid[ba$var < qf[q]]}
  names(lci) <- paste0("ci:", names(qf))
  lmvp[[names(lfilt)[t]]] <- lci
}
lmvp.abs <- lmvp
#-------------------------------------------
# Recurrent/cross-tissue low variance probes
#-------------------------------------------
# get recurrent probe ids from abs and bin methods
lmvp <- list()
for(t in names(lmvp.abs)){
  lii <- list(); li1 <- lmvp.abs[[t]]; li2 <- lmvp.bin[[t]]
  for(i in names(li1)){lii[[i]] <- unique(c(li1[[i]], li2[[i]]))}
  lmvp[[t]] <- lii
}
# get lowest recurrent abs filt probes
num.tissues <- 7
datlf <- lltx <- lmvp
txnames <- names(datlf)
txmvp <- matrix(nrow = 0, ncol = 2)
for(t in txnames){
  cgt <- lmvp.abs[[t]][[1]]
  dft <- data.frame(cgt, t, stringsAsFactors = F)
  txmvp <- rbind(txmvp, dft)
}
dft <- as.data.frame(table(txmvp[,1]))
dfft <- dft[dft[,2] == num.tissues,]
rt <- dfft
# get lowest recurrent bin filt probes
txmvp <- matrix(nrow = 0, ncol = 2)
for(t in txnames){
  cgt <- lmvp.bin[[t]][[1]]
  dft <- data.frame(cgt, t, stringsAsFactors = F)
  txmvp <- rbind(txmvp, dft)
}
dft <- as.data.frame(table(txmvp[,1]))
dfft <- dft[dft[,2] == 7,]
dim(dfft) # 553   2
# get unqiue low var probes
rt <- rbind(rt, dfft); rt <- rt[!duplicated(rt[,1]),] 
ptid <- rt[,1] # probe ids vector

#---------------------------------
# mean-filtered hypovar probe list
#---------------------------------
# get dnam summaries by probe of interest
meanrange.max <- 0.01
ldat.mean <- ldat.var <- list()
for(id in ptid){
  nr.mean <- nr.var <- c()
  for(t in names(lfilt)){
    lft <- lfilt[[t]]; lft <- lft[lft$cgid == id,]
    nr.mean <- c(nr.mean, lft$mean)
    nr.var <- c(nr.var, lft$var)
  }
  ldat.mean[[id]] <- nr.mean; ldat.var[[id]] <- nr.var
  # message(id)
}
did.mean <- do.call(rbind, ldat.mean); did.var <- do.call(rbind, ldat.var)
colnames(did.mean) <- colnames(did.var) <- names(lfilt)
for(c in seq(ncol(did.mean))){did.mean[,c] <- as.numeric(did.mean[,c])}
for(c in seq(ncol(did.var))){did.var[,c] <- as.numeric(did.var[,c])}
# apply dnam means filter
mean.min <- as.numeric(apply(did.mean, 1, min))
mean.max <- as.numeric(apply(did.mean, 1, max))
mean.range <- mean.max - mean.min
cgid.filt <- which(mean.range < meanrange.max)
ptidf <- ptid[cgid.filt]
save(ptidf, file = "cgids_lowvar-lowmean_nct7.rda")
