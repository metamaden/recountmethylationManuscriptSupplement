#!/usr/bin/env R

# Author: Sean Maden
# This script shows how to generate the main quality signal metrics (Table S2). 
# Note this requires a large file download to obtain the RGChannelSet object. 
# From this object are derived the, BeadArray control signals using `bactrl()`, 
# log2 medians of methylated and unmethylated signals, and the genotype 
# predictions and replicate inferences using `call_genotypes()` and 
# `check_snp_agreement()` from ewastools.

library(recountmethylationManuscriptSupplement)
library(recountmethylation)
library(HDF5Array)
library(minfi)
library(ewastools)

#----------
# load data
#----------
# download and load the H5SE red/grn signals dataset
# rgset <- recountmethylation::getdb_h5se_rg()
# get metadata
pkgname <- "recountmethylationManuscriptSupplement"
tables.dir <- system.file("extdata", "tables", package = pkgname)
tname <- "table-s1_mdpost_all-gsm-md.csv"
md <- read.csv(file.path(tables.dir, tname), header = TRUE)

#---------------------
# get quality controls
#---------------------
# beadarray control signals
rg <- getRed(rgset)
gg <- getGreen(rgset)
cdf <- as.data.frame(getProbeInfo(rgset, type = "Control"))
basignals <- bactrl(rs = t(rg), gs = t(gg), cdf = cdf)

# get log2 medians m/u signals
mset <- preprocessRaw(rgset)
ms <- getMeth(mset)
us <- getUnmeth(mset)
meth.l2med <- apply(ms, 2, function(x){log2(median(x))})
unmeth.l2med <- apply(us, 2, function(x){log2(median(x))})

# get pvalues and cutoffs
pvals <- minfi::detectionP(rgset)
detp.001 <- apply(pvals, 2, function(x){length(x[x > 0.01])})
detp.05 <- apply(pvals, 2, function(x){length(x[x > 0.05])})
detp.01 <- apply(pvals, 2, function(x){length(x[x > 0.1])})

#-------------------------------------
# detecting replicates from genotypes
#-------------------------------------
# get snp dnam
snp1.info <- getProbeInfo(rgset, type = "SnpI")
snp2.info <- getProbeInfo(rgset, type = "SnpII")
snp.addr <- c(snp1.info$AddressA, snp1.info$AddressB,
              snp2.info$AddressA)
beta.snp <- getSnpBeta(rgset)
# determine likely replicates by study id
athresh <- 0.1 # minimum probability
for(id in unique(md$gseid)){
 gsmv <- md[md$gseid == id,]$gsm
 sbi <- beta.snp[, colnames(beta.snp) %in% gsmv, drop = F]
 sgenoi <- call_genotypes(sbi, learn = F)
 sagree <- check_snp_agreement(sgenoi, donor_ids = colnames(sbi), 
                               sample_ids = colnames(sbi))
 ama = matrix(nrow = 0, ncol = 3)
 colnames(ama) <- c("gsm", "gseid", "cgsnp.gsm.geno")
 for(sai in sagree){
   am = data.frame(gsm = unique(c(sai$sample1, sai$sample2)), 
                   stringsAsFactors = F)
   am$gseid <- id; am$gsm.geno <- "NA"
   for(gsm in am$gsm){
     sai.cond <- (sai$sample1==gsm | sai$sample2==gsm) & 
       sai$agreement > athresh
     sai.gsm <- sai[sai.cond,]
     sai.id <- unique(c(sai.gsm$sample1, sai.gsm$sample2))
     am[am$gsm == gsm,]$gsm.geno <- paste(sai.id, collapse = ";")
   }
   ama <- rbind(ama, am)
 }
}
# append number of replicates
ama$num.shared <- unlist(lapply(ama$gsm.geno, function(x){
  length(unique(unlist(strsplit(x, ";"))))
}))

#---------------
# make new table
#---------------
# ba signals
mdf <- md[md$gsm %in% basignals$gsm,]
mdf <- mdf[order(match(mdf$gsm, basignals$gsm)),]
identical(mdf$gsm, basignals$gsm)
qcmd <- data.frame(gsm = basignals$gsm, gseid = mdf$gseid)
qcmd <- cbind(qcmd, basignals[,c(2:18)])
colnames(qcmd)[3:19] <- paste0("ba.", colnames(qcmd)[3:19])
# mu log2 medians
meth.l2med <- meth.l2med[order(match(names(meth.l2med), qcmd$gsm))]
unmeth.l2med <- unmeth.l2med[order(match(names(unmeth.l2med), qcmd$gsm))]
qcmd$meth.l2med <- meth.l2med
qcmd$unmeth.l2med <- unmeth.l2med
# replicate calls
ama <- ama[order(match(ama$gsm, qcmd$gsm)),]
qcmd$cgsnp.gsm.geno <- ama$gsm.geno
qcmd$cgsnp.nshared <- ama$num.shared

#----------------
# write new table
#----------------
#tname <- "table-s2_qcmd-allgsm.csv"
#write.table(tname, sep = ",")
