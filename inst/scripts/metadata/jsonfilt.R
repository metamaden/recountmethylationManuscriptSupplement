library(readr)
library(jsonlite)

# preamable, setting up shared variables
syssep = "/"
readpath = "gsm_json"
destpath = "gsm_json_filt"
filestem = ".json.filt"

lf.json <- list.files("gsm_json")
keys.list <- c("!Sample_characteristics_ch1",
               "!Sample_source_name_ch1",
               "!Sample_title")

# read json data, write filtered json data
for(i in 1:length(lf.json)){
  fni <- lf.json[i]
  gsmi <- unlist(strsplit(fni,"\\."))[2]
  writefn <- paste0(gsmi, filestem, collapse=".")
  writepath <- paste0(destpath, syssep, writefn,collapse="")
  rjsoni <- jsonlite::fromJSON(paste0(readpath,syssep,fni)) # read in json file
  
  # filter on valid sample-specific keys
  message("filtering keys for file ",i)
  rjsoni.keys <- colnames(rjsoni)
  rf <- list()
  for(k in 1:length(keys.list)){
    rekf <- rjsoni.keys[grepl(keys.list[k],rjsoni.keys)]
    for(f in 1:length(rekf)){
      rf[[rekf[f]]] <- as.character(unlist(rjsoni[rekf[f]]))
    }
  }
  
  # write properly formatted json with OUTSIDE BRACKETS (required!)
  message("writing filtered json data for file ",i)
  jsoni <- jsonlite::toJSON(rf, pretty=T, auto_unbox = T)
  write_lines("[", writepath)
  write_lines(jsoni, writepath, append=T)
  write_lines("]", writepath, append=T)
  
  message("finished file ",i)
}