#!user/bin/env R

# Author: Sean Maden
# This script defines a function, `soft2json`, that converts GSE 
# SOFT files into separate JSON-formatted samples/GSMs metadata files. 
# This function may be called from command line, eg:
# 
# > R soft2json.R <SOFT_FILENAME> <DESTINATION_DIRNAME>
# 
# JSON metadata files are used downstream for running the 
# MetaSRA-pipeline, after applying a filter for GSE-level information 
# (see the script `jsonfilt.R`)
#
# Arguments to the function `soft2json` are as follows:
# * gsmsoft_fn_list : List of GSE SOFT filenames to read in (vector).
# * gsm_json_destdir : Directory name to save new GSM JSON files (character).

require(jsonlite)

soft2json <- function(gsmsoft_fn_list, gsm_json_destdir){
  for(fn in gsmsoft_fn_list){
    gsmsoft_fn <- fn; linefile <- read.table(gsmsoft_fn,sep="\n")
    dffile <- as.data.frame(matrix(nrow=1,ncol=nrow(linefile)))
    colnames(dffile) <- gsub(" =.*","",linefile[,1])
    dffile[1,] <- gsub("!.*= ","",linefile[,1])
    jsoni <- toJSON(dffile, pretty=T); new.fn <- paste0(gsmsoft_fn,".json")
    new.fpath <- file.path(gsm_json_destdir,new.fn)
    write(jsoni,file = new.fpath)
  }
}

soft2json(gsmsoft_fn_list = commandArgs(T)[1], gsm_json_destdir = commandArgs(T)[2])


