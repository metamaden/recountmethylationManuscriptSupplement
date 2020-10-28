#!/usr/bin/env R

# Metadata term/label enrichment by storage condition quality 
# performance (BeadArray controls & M/U log2 median signals).

library(data.table)
library(recountmethylationManuscriptSupplement)

#----------
# load data
#----------
pkgname <- "recountmethylationManuscriptSupplement"
tables.dir <- system.file("extdata", "tables", package = pkgname)
# load qc data
fn1 <- "table-s2_qcmd-allgsm.csv"
qcmd <- fread(file.path(tables.dir, fn1), sep = ",", data.table = FALSE)
# load metadata
fn2 <- "table-s1_mdpost_all-gsm-md.csv"
md <- fread(file.path(tables.dir, fn2), sep = ",", data.table = FALSE)

#---------------------
# get filtered samples
#---------------------
mdf <- md[!md$storage == "NA",]
gsmv <- mdf$gsm
qcf <- qcmd[qcmd$gsm %in% gsmv,]
bas <- qcf[,c(1, 3:19)]
colnames(bas)[2:18] <- gsub("^ba\\.", "", colnames(bas)[2:18])
bat <- bathresh(bas)

#----------------
# helper function
#----------------
get_btm <- function(bgv, testv){
  bm <- matrix(nrow = 0, ncol = 3) # return matrix
  testv <- testv[!testv == "NA"] # na filt
  for(t in unique(testv)){
    pbg <- length(bgv[bgv == t])/length(bgv)
    pte <- length(testv[testv == t])/length(testv)
    pval <- binom.test(x = length(testv[testv == t]),
                       n = length(testv), p = pbg)$p.value
    et <- ifelse(pte > pbg, "overenriched", "underenriched")
    bti <- matrix(c(t, pval, et), nrow = 1)
    bm <- rbind(bm, bti)
  }
  colnames(bm) <- c("label", "pval", "enrichment.type")
  bm <- as.data.frame(bm, stringsAsFactors = FALSE)
  bm$pbh <- p.adjust(bm$pval, method = "BH")
  return(bm)
}

#----------------------
# m/u signal enrichment
#----------------------
# ffpe
gsm.ffpe <- mdf[grepl("FFPE", mdf$storage),]$gsm
ms1 <- mdf[mdf$gsm %in% gsm.ffpe,]
qs1 <- qcf[qcf[,1] %in% gsm.ffpe,c(1, 20, 21)]
sig.cond <- qs1[,2] < 11 & qs1[,3] < 11
msf1 <- mdf[mdf$gsm %in% qs1[sig.cond,]$gsm,]
stm.ffpe <- get_btm(bgv = unlist(strsplit(ms1$tissue, ";")),
                    testv = unlist(strsplit(msf1$tissue, ";")))
sdm.ffpe <- get_btm(bgv = unlist(strsplit(ms1$disease, ";")),
                    testv = unlist(strsplit(msf1$disease, ";")))

# ff
gsm.ff <- mdf[grepl("frozen", mdf$storage),]$gsm
ms2 <- mdf[mdf$gsm %in% gsm.ff,]
qs2 <- qcf[qcf[,1] %in% gsm.ff,c(1, 20, 21)]
sig.cond <- qs2[,2] < 11 & qs2[,3] < 11
msf2 <- mdf[mdf$gsm %in% qs2[sig.cond,]$gsm,]
stm.ff <- get_btm(bgv = unlist(strsplit(ms2$tissue, ";")),
                  testv = unlist(strsplit(msf2$tissue, ";")))
sdm.ff <- get_btm(bgv = unlist(strsplit(ms2$disease, ";")),
                  testv = unlist(strsplit(msf2$disease, ";")))

# overenriched labels -- tissue
pbh.max <- 1e-3
filt.cond <- stm.ffpe$enrichment.type == "overenriched" & stm.ffpe$pbh < pbh.max
stm.ffpe[filt.cond, 1]
# [1] "colorectal" "intestine"  "colon"      "rectum"     "mucosa"     "tumor"
filt.cond <- stm.ff$enrichment.type == "overenriched" & stm.ff$pbh < pbh.max
stm.ff[filt.cond, 1]
# [1] "nasal"            "epithelial"       "endocrine_system" "pancreas"

# overenriched labels -- disease
filt.cond <- sdm.ffpe$enrichment.type == "overenriched" & sdm.ffpe$pbh < pbh.max
sdm.ffpe[filt.cond, 1]
# [1] "case"    "normal"  "healthy" "control"
filt.cond <- sdm.ff$enrichment.type == "overenriched" & sdm.ff$pbh < pbh.max
sdm.ff[filt.cond, 1]
# character(0)

#----------------------------
# beadarray, label enrichment
#----------------------------
# FFPE samples
b1 <- bat[bat[,1] %in% mdf[grepl("FFPE", mdf$storage),]$gsm,]
m1 <- mdf[mdf$gsm %in% b1$gsm,]
min.failed <- 1
num.failed <- apply(b1, 1, function(x){length(x[x=="FAIL"])})
m1f <- m1[num.failed >= min.failed,]
btm.ffpe <- get_btm(bgv = unlist(strsplit(m1$tissue, ";")),
               testv = unlist(strsplit(m1f$tissue, ";")))
bdm.ffpe <- get_btm(bgv = unlist(strsplit(m1$disease, ";")),
               testv = unlist(strsplit(m1f$disease, ";")))
# FFPE, significantly enriched labels
pbh.max <- 1e-3
bcond <- btm$enrichment.type == "overenriched" & btm$pbh < pbh.max
btm[bcond, 1]
# [1] "metastasis"         "respiratory_system" "lung"               "colorectal"        
# [5] "intestine"          "colon"              "rectum"             "mucosa"            
# [9] "breast"
bcond <- bdm$enrichment.type == "overenriched" & bdm$pbh < pbh.max
bdm[bcond, 1]
# [1] "normal"        "healthy"       "control"       "lung_cancer"   "case"         
# [6] "breast_cancer"

# FF samples
b2 <- bat[bat[,1] %in% mdf[grepl("frozen", mdf$storage),]$gsm,]
m2 <- mdf[mdf$gsm %in% b2$gsm,]
min.failed <- 1
num.failed <- apply(b2, 1, function(x){length(x[x=="FAIL"])})
m2f <- m2[num.failed >= min.failed,]
btm.ff <- get_btm(bgv = unlist(strsplit(m2$tissue, ";")),
                    testv = unlist(strsplit(m2f$tissue, ";")))
bdm.ff <- get_btm(bgv = unlist(strsplit(m2$disease, ";")),
                    testv = unlist(strsplit(m2f$disease, ";")))
# FF, significantly enriched labels
pbh.max <- 1e-3
pfilt <- btm.ff$pbh < pbh.max & btm.ff$enrichment.type == "overenriched"
btm.ff[pfilt, 1]
# [1] "whole_blood" "sperm"       "oral"
pfilt <- bdm.ff$enrichment.type == "overenriched" & bdm.ff$pbh < pbh.max
bdm.ff[pfilt, 1]
# [1] "normal"  "control"
