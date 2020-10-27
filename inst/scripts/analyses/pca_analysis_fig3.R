#!/usr/bin/env R

# Author: Sean Maden
# Perform principal component anlayses (PCA) of autosomal DNAm
# across all samples, and within sample subsets. Inputs to this
# script are the large feature-hashed tables of each respective 
# DNAm dataset (access these from 
# https://recount.bio/data/recountmethylation_manuscript_supplement/data/).

library(data.table)
library(recountmethylation)
library(ggplot2)

#----------
# load data
#----------
pkgname <- "recountmethylationManuscriptSupplement"
pca.dir <- system.file("extdata", "pcadata", package = pkgname)
nct7.dir <- system.file("extdata", "nct7", package = pkgname)

# results data filenames
pca.all.fn <- "pcadat_allsamples.rda"; pca.noblood.fn <- "pcadat_all-noblood.rda"
pca.nct7.fn <- "pcadat_nct7.rda"; pca.nct6.fn <- "pcadat_nct6-nosperm.rda"

# load feature hashed data, and match to md
fhdat <- fread("fh1000_noobbeta-autosome.table", header = TRUE, sep = " ", data.table = FALSE)
fhdat <- fhdat[!duplicated(fhdat[,1]),]; rownames(fhdat) <- gsub("\\..*", "", fhdat[,1])

# get metadata with labels
mdfn <- "metadata_pca-labels.rda"
md <- get(load(file.path(pca.dir, mdfn)))

# nct 7 tissues
fn <- "gsmid-nct7.rda"
difid.nct7 <- get(load(file.path(nct7.dir, fn)))

#----------------
# format metadata
#----------------
# match md samples to feature hashed dataset
fhdat <- fhdat[,c(2:ncol(fhdat))]; gsm.int <- intersect(md$gsm, rownames(fhdat))
fhdat <- fhdat[rownames(fhdat) %in% gsm.int,]
fhdat <- fhdat[order(match(rownames(fhdat), gsm.int)),]
mdf <- md[md$gsm %in% gsm.int,]; mdf <- mdf[order(match(mdf$gsm, rownames(fhdat))),]
identical(rownames(fhdat), mdf$gsm)

#------------------------------------
# run pcas, save results data objects
#------------------------------------
set.seed(2)

# all samples
pca.all <- prcomp(fhdat)
pca.all$x <- pca.all$x[,c(1:2)]
save(pca.all, file = pca.all.fn)

# all, no blood/leukemia
gsm.filt <- md[!md$sample.label %in% c("blood", "leukemia"),]$gsm
fhdat.filt <- fhdat[rownames(fhdat) %in% gsm.filt,]; dim(fhdat.filt)
pca.noblood <- prcomp(fhdat.filt)
pca.noblood$x <- pca.noblood$x[,c(1:2)]
save(pca.noblood, file = pca.noblood.fn)

# nct7 with sperm
which.samp <- rownames(fhdat) %in% difid.nct7$gsmid
fhdat.filt <- fhdat[which.samp,]
pca.nct7 <- prcomp(fhdat.filt); save(pca.nct7, file = pca.nct7.fn)

# nct6 no sperm
dfid.nct6 <- difid.nct7[!difid.nct7$tissue == "sperm",]
which.samp <- rownames(fhdat) %in% dfid.nct6$gsmid
fhdat.filt <- fhdat[which.samp,]
pca.nct6 <- prcomp(fhdat.filt); save(pca.nct6, file = pca.nct6.fn)
