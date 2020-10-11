#!/usr/bin/env R

# Author: Sean Maden
# High variance autosomal DNAm analysis in 7 noncancer tissues/ This script 
# describes variance analysis on the ANOVA-filtered probe data from 7 tissues. 
# This takes as input the list object `lfilt` containing summary statistics for
# filtered cg probes in each tissue.

library(minfiData)

#----------
# load data
#----------
# get minfi-formatted probe annotations
cga <- getAnnotation(get(data("MsetEx")))
# get filtered cg datasets by tissue
pkgname <- "recountmethylationManuscriptSupplement"
nct7.dir <- system.file("extdata", "nct7", package = pkgname)
lf.fn <- "lfilt_cgtables_nct.rda"
lfilt <- get(load(file.path(nct7.dir, lf.fn)))

#---------------------------------------
# use bin method to get probes by tissue
#---------------------------------------
# params
qiv <- seq(0, 1, 0.01) # quantile filter
qwhich <- c(81, 86, 91, 96, 100) # c(17, 18, 19, 20)
binv <- seq(0, 1, 0.01)[1:100] # binned bval mean
# iter on ncts
lmvp <- list()
for(t in 1:length(lfilt)){
  ba <- lfilt[[t]]
  lci <- list(c(), c(), c(), c(), c())
  # iterate on betaval bins
  for(b in binv){
    bf <- ba[ba$mean >= b & ba$mean < b + 0.01, ] # filtered on mean bin
    qf <- quantile(bf$var, qiv)[qwhich]
    # iterate on ci values
    for(q in 1:length(qf)){
      lci[[q]] <- c(lci[[q]], bf$cgid[bf$var > qf[q]])
    }
    message(t, ":", b)
  }
  names(lci) <- paste0("ci:", names(qf))
  lmvp[[names(lfilt)[t]]] = lci
}
lmvp.bin <- lmvp

#---------------------------------------------------
# use absolute cutoff method to get probes by tissue
#---------------------------------------------------
lmvp <- list()
qiv <- seq(0, 1, 0.01) # quantile filter
qwhich <- c(81, 86, 91, 96, 100) # c(17, 18, 19, 20)
for(t in 1:length(lfilt)){
  ba <- lfilt[[t]]
  lci <- list(c(), c(), c(), c(), c())
  qf <- quantile(ba$var, qiv)[qwhich]
  for(q in 1:length(qf)){
    lci[[q]] <- ba$cgid[ba$var > qf[q]]
  }
  names(lci) <- paste0("ci:", names(qf))
  lmvp[[names(lfilt)[t]]] <- lci
}
lmvp.abs <- lmvp

#----------------------------------------------
# get probes with high tissue-specific variance
#----------------------------------------------
# collapse the lmvp lists
ptmax <- 7 # num shared tx for pan-tissue cgid
ncgmax <- 1000 # max probes by methods
lmvp <- list()
for(t in names(lmvp.abs)){
  lii <- list(); li1 <- lmvp.abs[[t]]; li2 <- lmvp.bin[[t]]
  for(i in names(li1)){lii[[i]] <- unique(c(li1[[i]], li2[[i]]))}
  lmvp[[t]] <- lii
}
datlf <- lltx <- lmvp; txnames <- names(datlf)
txmvpl <- list(); fnstem <- "filt"
bpdf1 <- bpdf5 <- bpdf10 <- matrix(nrow = 0, ncol = 3)
ltx <- c(); ltx2 <- c(); mvpcutl <- list()
for(wi in c(5, 4, 3)){
  which.lltx <- wi; txmvp <- matrix(nrow = 0, ncol = 2)
  for(tx in txnames){
    mvp.id <- lltx[[tx]][[which.lltx]]
    txmvp <- rbind(txmvp, data.frame(mvp.id, tx, stringsAsFactors = F))
  }
  # get pan-tissue probes
  dff <- as.data.frame(table(txmvp$mvp.id));
  ptid <- dff[dff[,2] == ptmax, 1]; num.ptid <- length(ptid)
  # filt for tx-specific probes only
  nontxid <- c(txmvp[,1][duplicated(txmvp[,1])]) # cgids in >1 tx
  txid <- txmvp[,1][!txmvp[,1] %in% nontxid]
  txmvp <- txmvp[txmvp[,1] %in% txid,]
  # bp data
  for(tx in txnames){
    mvpidt <- lltx[[tx]][[which.lltx]]
    mvpcutl[[tx]] <- mvpidt[mvpidt %in% txid]
    datt <- table(mvpidt %in% txid)
    nm <- matrix(c((as.numeric(datt)[1] - num.ptid), "Non-specific, other", tx,
                  as.numeric(datt)[2], "Tissue-specific", tx,
                  as.numeric(num.ptid), "Pan-tissue", tx), 
                nrow = 3, byrow = T)
    if(which.lltx == 5){bpdf1 <- rbind(bpdf1, nm)}
    if(which.lltx == 4){bpdf5 <- rbind(bpdf5, nm)}
    if(which.lltx == 3){bpdf10 <- rbind(bpdf10, nm)}
  }
  txmvp <- as.data.frame(txmvp, stringsAsFactors = F)
  colnames(txmvp) <- c("cgmvp", "txname")
  txmvpl[[length(txmvpl) + 1]] <- txmvp
}
names(txmvpl) = c("top1.above99.var", "top5.above95.var", "top10.above90.var")
txmvpl <- list(); fnstem <- "filt"
bpdf1 <- bpdf5 <- bpdf10 <- matrix(nrow = 0, ncol = 3)
ltx <- c(); ltx2 <- c(); mvpcutl <- list()
wi <- 1; txmvp <- matrix(nrow = 0, ncol = 2)
# get all mvp cgids
for(tx in txnames){
  mvp.id <- lltx[[tx]][[wi]]
  mvpdf <- data.frame(mvp.id, tx, stringsAsFactors = F)
  txmvp <- rbind(txmvp, mvpdf)
}
# get the tx specific cgids
dff <- as.data.frame(table(txmvp$mvp.id))
cgid.tx <- dff$Var1[dff$Freq == 1]
# tx top 1k var cgids, 2 methods
txfiltl <- list(); mt <- txmvpl[[2]]
for(tx in txnames){
  ltx <- lfilt[[tx]]
  # abs method filt
  txmvp.all <- lmvp.abs[[tx]][[1]]
  tx.mvpid <- txmvp.all[txmvp.all %in% cgid.tx]
  ltxf <- ltx[ltx$cgid %in% tx.mvpid,]
  abs.mvp.cgid <- ltxf[rev(order(ltxf$var)),]$cgid[1:ncgmax]
  # bin method filt, excluding abs method cgids
  txmvp.all <- lmvp.bin[[tx]][[1]]
  tx.mvpid <- txmvp.all[txmvp.all %in% cgid.tx & !txmvp.all %in% abs.mvp.cgid]
  ltxf <- ltx[ltx$cgid %in% tx.mvpid,]
  bin.mvp.cgid <- ltxf[rev(order(ltxf$var)),]$cgid[1:ncgmax]
  # append cgids
  txfiltl[[tx]] <- ltx[ltx$cgid %in% c(bin.mvp.cgid, abs.mvp.cgid),]
  message(tx)
}

ptid <- do.call(rbind, lapply(names(txfiltl), function(x){
  data.frame(tissue = rep(x, 2000), cgid = txfiltl[[x]][,1], stringsAsFactors = FALSE)
}))

#-------------
# save results
#-------------
txfiltl.fn <- "ltxfilt-dnamsummary_highvar_nct7.rda"
id.fn <- "cgids_highvar_nct7.rda"

save(txfiltl, file = txfiltl.fn)
save(ptid, file = id.fn)


