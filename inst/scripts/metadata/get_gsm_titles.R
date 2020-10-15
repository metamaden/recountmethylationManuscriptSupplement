#!/usr/bin/env R

# Author: Sean Maden
# This script describes how sample/GSM titles are extracted from the
# JSON-formatted metadata files. 

require(rjson)

#-------------
# specify path
#-------------
sfp <- "gsm_json_filt"
lf <- list.files(sfp)
lff <- lf[grepl(".*\\.json\\.filt$",lf)]

#------------------
# get sample titles
#------------------
gsmtitledf <- data.frame(matrix(c("", ""), nrow=0, ncol=2))
for(j in 1:length(lff)){
  gsm.id <- gsub("\\..*","",lff[j])
  json.lines <- readLines(paste0(sfp, "/", lff[j]))
  json.convert <- fromJSON(paste(json.lines,collapse="")) 
  st.catch <- as.character(unlist(json.convert)["!Sample_title"])
  stm <- matrix(c(gsm.id, st.catch), nrow=1, ncol=2)
  df.gsm <- as.data.frame(stm, stringsAsFactors = FALSE)
  gsmtitledf <- rbind(gsmtitledf, df.gsm)
  message(j)
}
colnames(gsmtitledf) <- c("gsm", "gsm_title")

#----------------
# save title data
#----------------
new.fn <- "gsm_jsontitledf.rda"
save(gsmtitledf, file = new.fn)
