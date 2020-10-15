#!/usr/bin/env R

# Author: Sean Maden
# This script describes the metadata preprocessing. This includes 
# harmonization of the mined metadata by moving terms under several 
# common variables for sample type, disease state, gender, age, anatomic 
# location, etc. It takes as input a list of available GSE SOFT file 
# data as a flat file. For each study, available data is assigned to 
# one of several variables.
#
# As input, this script loads the list of study-wise SOFT metadata 
# tables downloaded from GEO (see the script "make_gse_annolist.R" for 
# details). For info about acquiring the SOFT metadata, see the repo
# "recount-methylation-server".

#----------
# load data
#----------
# designate gse_tables path
dir <- "gse_tables"; dpath <- file.path(dir)
# load a list of the soft-derived terms
pkgname <- "recountmethylationManuscriptSupplement"
md.dir <- system.file("extdata", "metadata", package = pkgname) 
ldat.fn <- "geo_gse-atables_list.rda"
tls <- get(load(file.path(md.dir, ldat.fn)))
# load gsm titles
titles.fn <- "gsm_jsontitledf.rda"
tdf <- get(load(file.path(md.dir, titles.fn)))

#-------------------------
# process study/gse tables
#-------------------------
# begin main df
gat.all <- matrix(nrow = 0, ncol = 8)
colnames(gat.all) <- c("gsm","gse","sample_type","disease_state",
                       "gender","age","anatomic_location","misc")
# append tables data
fnv <- list.files(dpath)
for(gsei in seq(length(fnv))){
  source(file.path(dir, fnv[gsei]))
  gat.all <- rbind(gat.all, gat)
  # message(gsei)
}
gat.all <- gat.all[!duplicated(gat.all[,1]),]

#------------------
# append gsm titles
#------------------
# check for outersect gsm ids
d1 <- gat.all; d2 <- gsmtitledf
gsm.all <- unique(c(d1[,1], d2[,1]))
gsm1 <- gsm.all[!gsm.all %in% d1[,1]]
gsm2 <- gsm.all[!gsm.all %in% d2[,1]]
# append na slices as necessary
if(length(gsm1) > 0){
  nav <- rep(rep("NA", length(gsm1)), ncol(d1) - 1)
  mna <- matrix(c(gsm1, nav), nrow = length(gsm1), ncol = ncol(d1))
  d1 <- rbind(d1, mna)
}
if(length(gsm2) > 0){
  nav <- rep(rep("NA", length(gsm2)), ncol(d2) - 1)
  mna <- matrix(c(gsm2, nav), nrow = length(gsm2), ncol = ncol(d2))
  d2 <- rbind(d2, mna)
}
# reorder and assign title var
match.gsm1 <- match(as.character(d1[,1]), as.character(d2[,1]))
order.gsm1 <- order(match.gsm1)
d1 <- d1[order.gsm1,]
match.gsm2 <- match(as.character(d2[,1]), as.character(d1[,1]))
order.gsm2 <- order(match.gsm2)
d2 <- d2[order.gsm2,]
cond <- identical(as.character(d2[,1]), as.character(d1[,1]))
if(cond){
  d1 <- as.data.frame(d1, stringsAsFactors = FALSE)
  d1$gsm_title <- as.character(d2[,2])
}
mdpre <- d1

#--------------------------------
# save the postprocessed metadata
#--------------------------------
mdpre.fn <- "md_preprocess"
save(mdpre, file = paste0(mdpre.fn, ".rda"))
write.csv(mdpre, file = paste0(mdpre.fn, ".csv"))
