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
# load a list of the soft-derived terms
pkgname <- "recountmethylationManuscriptSupplement"
md.dir <- system.file("extdata", "metadata", package = pkgname) 
ldat.fn <- "geo_gse-atables_list.rda"
tls <- tls <- get(load(file.path(md.dir, ldat.fn)))

#-------------------------
# process study/gse tables
#-------------------------
# begin main df
gat.all <- matrix(nrow = 0, ncol = 8)
colnames(gat.all) <- c("gsm","gse","sample_type","disease_state",
                       "gender","age","anatomic_location","misc")
# append tables data
dir <- "gse_tables"
fnv <- list.files(dir)
for(gsei in seq(length(fnv))){
  source(file.path(dir, fnv[gsei]))
  gat.all <- rbind(gat.all, gat)
  message(gsei)
}

gat.all <- gat.all[!duplicated(gat.all[,1]),]

#--------------------------------
# save the postprocessed metadata
#--------------------------------
mdpre <- gat.all
mdpre <- as.data.frame(mdpre, stringsAsFactors = FALSE)

mdpre.fn <- "md_preprocess"
save(mdpre, file = paste0(mdpre.fn, ".rda"))
write.csv(mdpre, file = paste0(mdpre.fn, ".csv"))
