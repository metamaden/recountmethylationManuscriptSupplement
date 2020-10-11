#!/usr/bin/env R

# Author: Sean Maden
# Describes how to generate model based predictions for sample metadata, 
# from DNAm data, for age, sex, and blood cell type fractions. This script
# takes as input 2 database HDF5-SummarizedExperiment-type database files
# downloaded from the data server at https://recount.bio/data. After applying
# model-based predictions, the table `mdmod`, which is appended with other 
# mined and learned metadata in the final sample metadata table (Table S1).

library(recountmethylation)
library(wateRmelon)
library(minfi)

#----------
# load data
#----------
# download and load the H5SE red/grn signals dataset
rgset <- recountmethylation::getdb_h5se_rg()
gmset <- recountmethylation::getdb()

#----------------
# get predictions
#----------------
# sex predictions
sexpred <- getSex(gmset)
# cell type predictions
celltypepred <- estimateCellCounts(rgset)
# age predictions
grset <- preprocessNoob(rgset)
predage <- agep(getBeta(grset))

#-----------
# make table
#-----------
mdmod <- data.frame(gsm = rgset$gsm, stringsAsFactors = FALSE)
mdmod$predsex <- sexpred$predictedSex
mdmod$predage <- predage[,1]
mdmod <- cbind(mdmod, celltypepred)
