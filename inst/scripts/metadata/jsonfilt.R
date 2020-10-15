#!/usr/bin/env R

# Author: Sean Maden
# This script shows how JSON-formatted sample/GSM metadata, 
# extracted from GSE SOFT files, was filtered prior to running in
# MetaSRA-pipeline. Most likely/highest-confidence sample type 
# predictions from this run were stored under the "sampletype"
# column in the main metadata table (Table S1).

library(readr)
library(jsonlite)

#------------------------
# assign dir's and path's
#------------------------
# input dir/path
json.dn <- "gsm_json"
readpath <- file.path("recount-methylation-files", json.dn)

# output dir/path
jfilt.dn <- json.dn # paste0(json.dn, "_filt")
destpath <- file.path("recount-methylation-files", jfilt.dn)
if(!dir.exists(destpath)){dir.create(destpath)}

# new files stem
filestem <- "json.filt"

# keys of interest
keys.list <- c("!Sample_characteristics_ch1", "!Sample_source_name_ch1", 
               "!Sample_title")

#------------------
# filter json files
#------------------
# get unfiltered soft files
lf.json <- list.files(readpath)
lf.json <- lf.json[!grepl(".*filt.*", lf.json)]
# write filtered data as new json files
for(i in seq(length(lf.json))){
  fni <- lf.json[i]
  ts <- unlist(strsplit(fni,"\\."))[1] # timestamp
  gsmi <- unlist(strsplit(fni,"\\."))[2] # gsm id
  writefn <- paste(ts, gsmi, filestem, sep = ".")
  writepath <- file.path(destpath, writefn)
  rjsoni <- jsonlite::fromJSON(file.path(readpath,fni))
  # filter on valid sample-specific keys
  message("filtering keys for file ",i)
  rjsoni.keys <- colnames(rjsoni); rf <- list()
  for(k in keys.list){rekf <- rjsoni.keys[grepl(k, rjsoni.keys)]
    for(f in rekf){rf[[f]] <- as.character(unlist(rjsoni[f]))}
  }
  # write formatted json with top and bottom outside brackets
  message("writing filtered json data for file ",i)
  jsoni <- jsonlite::toJSON(rf, pretty=T, auto_unbox = T)
  write_lines("[", writepath); write_lines(jsoni, writepath, append=T)
  write_lines("]", writepath, append=T)
  message("finished file ",i)
}
