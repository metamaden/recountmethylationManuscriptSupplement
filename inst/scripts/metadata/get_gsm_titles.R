#!/usr/bin/env R

# Author: Sean Maden
# This script describes how sample/GSM titles are extracted from the
# JSON-formatted metadata files. 

require(rjson)

#-------------
# specify path
#-------------
dir <- "gsm_json"; dpath <- file.path("recount-methylation-files", dir)
lf <- list.files(dpath); ffilt <- grepl(".*\\.json\\.filt$", lf)
lff <- lf[ffilt]

#------------------
# get sample titles
#------------------
gsmtitledf <- matrix(nrow = 0, ncol = 2)
for(jf in lff){
  jsym <- unlist(strsplit(jf, "\\.")); gsm.id <- jsym[2]
  json.lines <- readLines(file.path(dpath, jf))
  json.convert <- fromJSON(paste(json.lines,collapse="")) 
  st.catch <- as.character(unlist(json.convert)["!Sample_title"])
  stm <- matrix(c(gsm.id, st.catch), nrow=1)
  gsmtitledf <- rbind(gsmtitledf, stm)
}
colnames(gsmtitledf) <- c("gsm", "gsm_title")

#----------------
# save title data
#----------------
new.fn <- "gsm_jsontitledf.rda"
save(gsmtitledf, file = new.fn)
