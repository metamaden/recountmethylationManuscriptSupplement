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

library(recountmethylationManuscriptSupplement)
library(HDF5Array)
library(recountmethylation)
library(wateRmelon) # for age predictions
library(minfi) # for cell type and sex predictions

#----------
# load data
#----------
# path to postprocessed metadata
pkgname <- "recountmethylationManuscriptSupplement"
md.dir <- system.file("extdata", "metadata", package = pkgname) 
md.fn <- "metadata_preprocessing.rda"
mdpre <- get(load(file.path(md.dir, md.fn)))

# download and load the H5SE red/grn signals dataset
#rgset <- recountmethylation::getdb_h5se_rg()
#grset <- recountmethylation::getdb_h5se_gr()

# load the delayed arrray databases
rgset.path <- "remethdb_h5se-rg_0-0-1_1590090412"
grset.path <- "remethdb-h5se_gr_0-0-1_1590090412"
rgset <- loadHDF5SummarizedExperiment(rgset.path)
grset <- loadHDF5SummarizedExperiment(grset.path)

#----------------
# get predictions
#----------------
# start mdmod
mdmod <- matrix(nrow = 0, ncol = 8)

# do pred on sample index blocks
blocks <- getblocks(ncol(rgset), 50) # get input indices
for(b in blocks){
  rgf <- rgset[, unlist(b)]
  rgf.matrix <- RGChannelSet(Green = as.matrix(getGreen(rgf)),
                             Red = as.matrix(getRed(rgf)),
                             annotation = annotation(rgf))
  celltypepred <- estimateCellCounts(rgf.matrix) # cell type predictions
  msf <- mapToGenome(preprocessRaw(rgf))
  sexpred <- getSex(msf) # sex predictions
  grf <- grset[,colnames(rgf)]
  predage <- agep(getBeta(grf)) # age predictions
  # append results
  mdf <- cbind(predage, cbind(sexpred[,3], celltypepred))
  mdmod <- rbind(mdmod, mdf)
  message("finished up to sample index ", max(unlist(b)))
}

# format table
mdmod <- as.data.frame(mdmod, stringsAsFactors = FALSE)
colnames(mdmod) <- c("predage", "predsex", 
                     paste0("predcell.", colnames(mdmod)[3:8]))
mdmod$gsm <- gsub("\\..*", "", rownames(mdmod))
mdmod <- mdmod[!duplicated(mdmod$gsm),]

#-----
# save
#-----
table.fn <- "mdmod"
save(mdmod, file = paste0(table.fn, ".rda"))
write.csv(mdmod, file = paste0(table.fn, ".csv"))
