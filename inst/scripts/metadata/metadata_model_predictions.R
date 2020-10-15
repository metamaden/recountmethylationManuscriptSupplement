#!/usr/bin/env R

# Author: Sean Maden
# Describes how to generate model based predictions for sample metadata, 
# from DNAm data, for age, sex, and blood cell type fractions. 
# 
# This script takes as input 2 database HDF5-SummarizedExperiment-type 
# database files downloaded from the data server at https://recount.bio/data. 
# After applying model-based predictions, the table `mdmod`, which is 
# appended with other mined and learned metadata in the final sample 
# metadata table (Table S1).

library(recountmethylation)
library(wateRmelon)
library(minfi)

#----------
# load data
#----------
# path to postprocessed metadata
pkgname <- "recountmethylationManuscriptSupplement"
md.dir <- system.file("extdata", "metadata", package = pkgname) 
md.fn <- "metadata_preprocessing.rda"
mdpre <- get(load(file.path(md.dir, md.fn)))

# download and load the H5SE red/grn signals dataset
rgset <- recountmethylation::getdb_h5se_rg()
gmset <- recountmethylation::getdb_h5se_gm()
grset <- recountmethylation::getdb_h5se_gr()

#----------------
# get predictions
#----------------
# sex predictions
sexpred <- getSex(gmset)
# cell type predictions
celltypepred <- estimateCellCounts(rgset)
# age predictions
predage <- agep(getBeta(grset))

#-----------
# make table
#-----------
mdmod <- data.frame(gsm = rownames(predage))
mdmod <- cbind(mdmod, cbind(sexpred[,3], cbind(predage, celltypepred)))
colnames(mdmod) <- c("gsm", "predsex", "predage", 
                     paste0("predcell.", colnames(mdmod)[4:9]))
