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
# specify index blocks for processes
blocks <- getblocks(slength = ncol(rgset), bsize = 50)
# get quality metrics for each index block
ms <- matrix(nrow = 0, ncol = 20)
cdf <- as.data.frame(getProbeInfo(rgset, type = "Control"))
for(b in blocks){
  rgf <- rgset[, b]; colnames(rgf) <- gsub("\\..*", "", colnames(rgf))
  # beadarray control signals
  redsignal <- getRed(rgf); greensignal <- getGreen(rgf)
  basignals <- bactrl(rs = t(redsignal), gs = t(greensignal), cdf = cdf)
  # get log2 medians m/u signals
  mset <- preprocessRaw(rgf)
  ms <- getMeth(mset); meth.l2med <- apply(ms, 2, function(x){log2(median(x))})
  us <- getUnmeth(mset); unmeth.l2med <- apply(us, 2, function(x){log2(median(x))})
  # bind qc signals
  mi <- cbind(basignals, 
              data.frame(meth.l2med = meth.l2med, 
                         unmeth.l2med = unmeth.l2med,
                         stringsAsFactors = FALSE))
  ms <- rbind(ms, mi)
}

#-------------------------------------
# detecting replicates from genotypes
#-------------------------------------
# get snp dnam
snp1.info <- getProbeInfo(rgset, type = "SnpI")
snp2.info <- getProbeInfo(rgset, type = "SnpII")
snp.addr <- c(snp1.info$AddressA, snp1.info$AddressB, snp2.info$AddressA)
beta.snp <- getSnpBeta(rgset)
colnames(beta.snp) <- gsub("\\..*", "", colnames(beta.snp))
# determine likely replicates by study id
athresh <- 0.1 # minimum probability
ama = matrix(nrow = 0, ncol = 3) # main identity matrix
colnames(ama) <- c("gsm", "gseid", "cgsnp.gsm.geno")
for(id in unique(md$gseid)){
 gsmv <- md[md$gseid == id,]$gsm
 sbi <- beta.snp[, colnames(beta.snp) %in% gsmv, drop = FALSE]
 sbi <- as.matrix(sbi); class(sbi) <- "numeric"
 sgenoi <- call_genotypes(sbi, learn = FALSE)
 sagree <- check_snp_agreement(sgenoi, donor_ids = colnames(sbi), 
                               sample_ids = colnames(sbi))
 for(sai in sagree){
   am <- data.frame(gsm = unique(c(sai$sample1, sai$sample2)), 
                   stringsAsFactors = FALSE)
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
 message("finished study ", id)
}
# append replicates count
ama$num.shared <- unlist(lapply(ama$gsm.geno, function(x){
  length(unique(unlist(strsplit(x, ";"))))
}))

# append replicates info
qcmd <- ms
ama <- ama[ama$gsm %in% qcmd$gsm,]
ama <- ama[order(match(ama$gsm, qcmd$gsm)),]
if(identical(ama$gsm, qcmd$gsm)){
  qcmd$cgsnp.gsm.geno <- ama$gsm.geno
  qcmd$cgsnp.nshared <- ama$num.shared
} else{stop("error matching gsm ids for qcmd, ama...")}

#----------------
# write new table
#----------------
tname <- "table-s2_qcmd-allgsm.csv"
write.csv(qcmd, file = tname)
