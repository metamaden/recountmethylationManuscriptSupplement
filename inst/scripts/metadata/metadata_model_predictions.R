#!/usr/bin/env R

# Describes how to generate model based predictions for 
# sample metadata, from DNAm data.

library(recountmethylation)
library(wateRmelon)
library(minfi)

#----------
# load data
#----------
# download and load the H5SE red/grn signals dataset
# rgset <- recountmethylation::getdb_h5se_rg()
# gmset <- recountmethylation::getdb()

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
