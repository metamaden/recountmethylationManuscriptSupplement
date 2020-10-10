#!/usr/bin/env R

# Describes how to generate model based predictions for 
# sample metadata, from DNAm data.

library(wateRmelon)
library(minfi)
library(FlowSorted.Blood.450k)
library(recountmethylation)
library(HDF5Array)

pkgname <- "recountmethylationManuscriptSupplement"
tables.dir <- system.file("extdata", "tables", package = pkgname)
tname <- "table-s1_mdpost_all-gsm-md.csv"
md <- read.csv(file.path(tables.dir, tname), header = TRUE)

#----------
# load data
#----------
# download and load the H5SE red/grn signals dataset
# rgset <- recountmethylation::getdb_h5se_rg()
# gmset <- recountmethylation::getdb()

rgpath <- "remethdb-h5_rg_1590090412_0-0-1.h5"
gsmv <- c("GSM1873723", "GSM2459564", "GSM3584264", "GSM2363655", "GSM1858832", "GSM2293090")
rgset <- getrg(dbn = rgpath, gsmv = gsmv)

gmset <- loadHDF5SummarizedExperiment("remethdb-h5se_gm_0-0-1_1590090412")

gsm.int <- intersect(md$gsm, gmset$gsm)
mdf <- md[md$gsm %in% gsm.int,]
gmset <- gmset[,gmset$gsm %in% gsm.int]
mdf <- mdf[order(match(mdf$gsm, gmset$gsm)),]
identical(gmset$gsm, mdf$gsm)

#----------------
# sex predictions
#----------------
sexpred <- getSex(gmset)
table(sexpred$predictedSex, gmset$predsex)
table(sexpred$predictedSex, mdf$predsex)

#----------------------
# cell type predictions
#----------------------
celltypepred <- estimateCellCounts(rgset)

mdf <- md[md$gsm %in% rgset$gsm,]
mdf <- mdf[order(match(mdf$gsm, rgset$gsm)),]
identical(mdf$gsm, rgset$gsm)

#----------------
# age predictions
#----------------

#----------
# save data
#----------