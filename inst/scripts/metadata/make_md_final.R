#!/usr/bin/env R

# Author: Sean Maden
# This script shows how the final metadata is assembled from various
# sources, including learned annotations/labels, sample type predictions
# from MetaSRA-pipeline, model-based predictions for age, sex, and blood
# cell types, and storage condition info. This corresponds to Table S1 
# from the main manuscript.

library(data.table)
library(recountmethylationManuscriptSupplement)

#----------
# load data
#----------
# path to metadata files
pkgname <- "recountmethylationManuscriptSupplement"
metadata.dir <- system.file("extdata", "metadata", package = pkgname)
# postprocessed metadata
mdfin <- mdpost <- get(load(file.path(metadata.dir, "md_postprocess.rda")))
# model based predictions
mdmod <- get(load(file.path(metadata.dir, "mdmod.rda")))
# storage procedure annotations
prepd <- get(load(file.path(metadata.dir, "prepmd.rda"))) 
# metasra-pipeline, mapped and predicted labels
mdmap <- get(load(file.path(metadata.dir, "mdmap-gsm_35k.rda"))) 
# formatted Cellosaurus records
# note, this file needs to be decompressed
ccf <- fread(file.path(metadata.dir, "ccformat.txt"), sep = " ", header = TRUE)
# rsheet object
rsh <- fread(file.path(metadata.dir, "1548876778.rsheet"), sep = ' ', 
             select = c(1:3, 6:8), header = TRUE, fill = TRUE)

#--------------------------------
# get sample type pred from mdmap
#--------------------------------
# append best type predictions
gsmid <- mdmap$gsmid
stype <- gsub("'", "",
              gsub(";.*", "",
                   gsub("^.*'sample type':", "", mdmap$msrap_flatjson)))
spred <- as.numeric(gsub("'", "",
                         gsub(";.*", "",
                              gsub("^.*'sample-type confidence':", "", 
                                   mdmap$msrap_flatjson))))
spred <- round(as.numeric(spred), digits = 3)
dfmd <- data.frame(gsm = gsmid, type = stype, pred = spred, stringsAsFactors = F)
dfmd <- dfmd[dfmd$gsm %in% mdfin$gsm,]
mdp.filt <- !mdfin$gsm %in% dfmd$gsm
na.val <- rep("NA", nrow(mdfin[mdp.filt,]))
dfp <- data.frame(gsm = mdfin[mdp.filt,]$gsm, type = na.val,
                  pred = na.val, stringsAsFactors = F)
dfmd <- rbind(dfmd, dfp)
dfmd <- dfmd[order(match(dfmd$gsm, mdfin$gsm)),]
mdfin$sampletype <- paste(paste("msraptype", dfmd$type, sep  =":"),
                          paste("msrapconf", dfmd$pred, sep  =":"),
                          sep = ";")

#---------------------------------------
# append cellosaurus info for cell lines
#---------------------------------------
mdst <- rep("NA", nrow(md)); miscdat <- md$misc
cll <- ifelse(grepl(".*cell_line.*", miscdat) & !grepl(".*cell_line:NA.*", miscdat),
              gsub("(;|$).*", "", gsub("(^|;).*cell_line: ", "", miscdat)),  "NA")
cll <- cll[!cll=="NA"]; cllf <- cll[cll %in% ccf$ID]
ccff <- ccf[ccf$ID %in% cllf,]
ccff$CA <- tolower(substr(ccff$CA, 4, nchar(ccff$CA)))
for(r in 1:nrow(ccff)){
  dati <- as.character(ccff[r,]); cf <- get_filt(get_pstr(dati[1]))
  if(length(cf[cf]) > 0){
    mdfin$sampletype <- appendvar("sampletype", "cell_line", cf)
    mdfin$sampletype <- appendvar("sampletype", paste0("ccid:", dati[1]), cf)
    mdfin$sampletype <- appendvar("sampletype", paste0("ccacc:", dati[2]), cf)
    mdfin$sampletype <- appendvar("sampletype", paste0("cccat:", dati[5]), cf)
  }
  message(r)
}

#-------------------------
# storage cond annotations
#-------------------------
# from prepd object
mdfin$storage <- "NA"
prepd$storage <- prepd$preparation
prepd$storage <- ifelse(!grepl(".*FFPE.*", prepd$storage),
                        "F;frozen", "FFPE;formalin_fixed_paraffin_embedded")
prepd <- prepd[,c("gsm", "storage")]; mdfin <- mdfin[,c("gsm", "storage")]
mdfin <- match1to2(as.matrix(mdfin), as.matrix(prepd), ci1 = 1, ci2 = 1)
mdfin <- as.data.frame(mdfin, stringsAsFactors = FALSE)

#----------------------
# append remaining data
#----------------------
# array files info
mdfin <- match1to2(as.matrix(mdfin), as.matrix(rsh), 1, 1)
# gseids from mdmap
mdfin <- match1to2(as.matrix(mdfin), as.matrix(mdmap[,c(1:2)]), 1, 1)
# match model predictions
mdfin <- match1to2(as.matrix(mdfin), as.matrix(mdmod), 1, 1)

mdfin <- as.data.frame(mdfin[,c(1:7, 10:13, 17:24)], stringsAsFactors = FALSE)

#-------------
# save objects
#-------------
# save data
mdfin <- "md_final"
save(mdfin, file = paste0(mdfin, ".rda"))

# save table
tname <- "table-s1_mdpost_all-gsm-md"
write.csv(mdfin, file = paste0(tname, ".csv"), row.names = FALSE)
