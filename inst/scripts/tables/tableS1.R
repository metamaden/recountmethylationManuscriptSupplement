#!/usr/bin/env R

# Author: Sean Maden
# This script describes how the sample metadata table is assembled from 
# various data files. GEO mined GSM titles, GSE IDs, and learned labels
# are contained in the `mdpost` object. MetaSRA-pipeline sample type 
# predictions for each GSM ID are included in the flat file `mdmap`. 
# Manually annotated storage conditions are contained in the object
# `prepd`. For details see `metadata_postprocessing.R` and other scripts
# for metadata preparations.

library(recountmethylationManuscriptSupplement)

#----------
# load data
#----------
pkgname <- "recountmethylationManuscriptSupplement"
metadata.dir <- system.file("extdata", "metadata", package = pkgname)
# postprocessed metadata
mdpost <- get(load(file.path(metadata.dir, "md-postprocess.rda")))
# storage procedure annotations
prepd <- get(load(file.path(metadata.dir, "prepmd.rda"))) 
# metasra-pipeline, mapped and predicted labels
mdmap <- get(load(file.path(metadata.dir, "mdmap-gsm_35k.rda"))) 
# formatted Cellosaurus records
# note, this file needs to be decompressed
ccf <- fread(file.path(metadata.dir, 'ccformat.txt'), sep = ' ', header = T) 

#-----------------------
# storage condition data
#-----------------------
# from prepd object
mdpost$storage <- "NA"
prepd$storage <- prepd$preparation
prepd$storage <- ifelse(!grepl(".*FFPE.*", prepd$storage),
                       "F;frozen", "FFPE;formalin_fixed_paraffin_embedded")
prepd <- prepd[,c("gsm", "storage")]; mf <- mdpost[,c("gsm", "storage")]
prepd <- rbind(prepd, mf[!mf$gsm %in% prepd$gsm,])
prepd <- prepd[order(match(prepd$gsm, mf$gsm)),]
mdpost$storage <- prepd$storage

#----------------------------------------------
# sample type predictions from metasra-pipeline
#----------------------------------------------
# append most likely sample type predictions
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
dfmd <- dfmd[dfmd$gsm %in% mdpost$gsm,]
mdp.filt <- !mdpost$gsm %in% dfmd$gsm
na.val <- rep("NA", nrow(mdpost[mdp.filt,]))
dfp <- data.frame(gsm = mdpost[mdp.filt,]$gsm, type = na.val,
                  pred = na.val, stringsAsFactors = F)
dfmd <- rbind(dfmd, dfp)
dfmd <- dfmd[order(match(dfmd$gsm, mdpost$gsm)),]
mdpost$sampletype <- paste(paste("msraptype", dfmd$type, sep  =":"),
                          paste("msrapconf", dfmd$pred, sep  =":"),
                          sep = ";")

# sample type, cellosaurus cell line info
mdst <- rep("NA", nrow(md)); miscdat <- md$misc
cll <- ifelse(grepl(".*cell_line.*", miscdat) & !grepl(".*cell_line:NA.*", miscdat),
             gsub("(;|$).*", "", gsub("(^|;).*cell_line: ", "", miscdat)),  "NA")
cll <- cll[!cll=="NA"]; cllf <- cll[cll %in% ccf$ID]
# filt cellosaurus
ccff <- ccf[ccf$ID %in% cllf,]
ccff$CA <- tolower(substr(ccff$CA, 4, nchar(ccff$CA)))
for(r in 1:nrow(ccff)){
  dati <- as.character(ccff[r,]); cf <- get_filt(get_pstr(dati[1]))
  if(length(cf[cf]) > 0){
    mdpost$sampletype <- appendvar("sampletype", "cell_line", cf)
    mdpost$sampletype <- appendvar("sampletype", paste0("ccid:", dati[1]), cf)
    mdpost$sampletype <- appendvar("sampletype", paste0("ccacc:", dati[2]), cf)
    mdpost$sampletype <- appendvar("sampletype", paste0("cccat:", dati[5]), cf)
  }
  message(r)
}

#-----------
# save table
#-----------
tname <- "table-s1_mdpost_all-gsm-md"
save(tname, file = paste0(tname, ".rda"))
write.csv(tname, file = paste0(tname, ".csv"), row.names = FALSE)
save(mdpost, file = nfn)
