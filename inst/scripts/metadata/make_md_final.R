#!/usr/bin/env R

# Author: Sean Maden
# This script shows how the final metadata is assembled from various
# sources, including learned annotations/labels, sample type predictions
# from MetaSRA-pipeline, model-based predictions for age, sex, and blood
# cell types, and storage condition info. This corresponds to Table S1 
# from the main manuscript.

library(recountmethylationManuscriptSupplement)

#----------
# load data
#----------
# path to metadata files
pkgname <- "recountmethylationManuscriptSupplement"
metadata.dir <- system.file("extdata", "metadata", package = pkgname)

# postprocessed metadata
mdpost <- get(load(file.path(metadata.dir, "md-postprocess.rda")))

# model based predictions
mdmod <- get(load(file.path(metadata.dir, "mdmod.rda")))

# storage procedure annotations
prepd <- get(load(file.path(metadata.dir, "prepmd.rda"))) 

# metasra-pipeline, mapped and predicted labels
mdmap <- get(load(file.path(metadata.dir, "mdmap-gsm_35k.rda"))) 

# formatted Cellosaurus records
# note, this file needs to be decompressed
ccf <- fread(file.path(metadata.dir, 'ccformat.txt'), sep = ' ', header = T)

#------------------
# gseids from mdmap
#------------------
mdfin <- match1to2()

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

# append cellosaurus info for cell lines
mdst <- rep("NA", nrow(md)); miscdat <- md$misc
cll <- ifelse(grepl(".*cell_line.*", miscdat) & !grepl(".*cell_line:NA.*", miscdat),
              gsub("(;|$).*", "", gsub("(^|;).*cell_line: ", "", miscdat)),  "NA")
cll <- cll[!cll=="NA"]; cllf <- cll[cll %in% ccf$ID]
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

#-------------------------
# storage cond annotations
#-------------------------
# from prepd object
mdfin$storage <- "NA"
prepd$storage <- prepd$preparation
prepd$storage <- ifelse(!grepl(".*FFPE.*", prepd$storage),
                        "F;frozen", "FFPE;formalin_fixed_paraffin_embedded")
prepd <- prepd[,c("gsm", "storage")]; mf <- mdfin[,c("gsm", "storage")]
prepd <- rbind(prepd, mf[!mf$gsm %in% prepd$gsm,])
prepd <- prepd[order(match(prepd$gsm, mf$gsm)),]
mdfin$storage <- prepd$storage

#------------------------
# match model predictions
#------------------------
mdfin <- match1to2()

#-----------
# save table
#-----------
tname <- "table-s1_mdpost_all-gsm-md"
save(tname, file = paste0(tname, ".rda"))
write.csv(tname, file = paste0(tname, ".csv"), row.names = FALSE)
save(mdpost, file = nfn)
