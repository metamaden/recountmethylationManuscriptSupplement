#!/usr/bin/env R

# Author: Sean Maden
#
# Make annotations list for GEO studies (GSE IDs). 
# This script reads in filtered JSON-formatted sample/GSM
# data, derived from GEO GSE SOFT files, and returns a
# list of mined metadata as a list. The input is a filtered
# table of IDs returned from a call to the script 
# "edirect_query.py", available here:
# https://github.com/metamaden/recount-methylation-server/blob/master/src/edirect_query.py

require(data.table)
require(rjson)

#---------------------------------
# load data, specify json location
#---------------------------------
# read in a query filt table 
eqf <- "gsequery_filt.1552256509"
x <- scan(eqf, what="", sep="\n")
gsel <- list()
for(i in 1:length(x)){
  ssi <- unlist(strsplit(x[[i]]," "))
  gsel[[ssi[1]]] <- ssi[2:length(ssi)]
  message(i)
}
# filtered json data location
dn <- "gsm_json_filt"
tgse.list <- list()

#---------------------
# make gse tables list
#---------------------
for(g in 1:length(gsel)){
  lgse <- list()
  gseid <- as.character(names(gsel)[g])
  ggse <- gsel[[g]]
  # parse json metadata
  for(j in 1:length(ggse)){
    fnj <- paste0("./", dn, "/", ggse[j], ".json.filt")
    if(file.exists(fnj)){
      lgse[[ggse[j]]] <- fromJSON(paste(readLines(fnj), collapse=""))
    }
    message(j)
  }
  if(length(lgse)>0){
    # make table columns
    tcols <- c()
    for(l in 1:length(lgse)){
      gsmid = names(lgse)[l]
      gsmval = unlist(lgse[[l]])
      for(k in 1:length(gsmval)){
        if(grepl(":",gsmval[k])){
          kk <- as.character(gsub(":.*","",gsmval[k]))
          if(!kk=="" & !kk %in% tcols){
            tcols <- c(tcols, kk) 
          }
        }
        #message(k)
      }
      message(l)
    }
    tcols <- c("gsm","gse",tcols)
    tgse <- matrix(nrow=0,ncol=length(tcols))
    colnames(tgse) <- tcols
    # coerce to study-specific table
    for(l in 1:length(lgse)){
      gsmid = names(lgse)[l]
      gsmval = unlist(lgse[[l]])
      gvk <- c(gsmid, gseid)
      # loop over gsm values
      for(c in 3:ncol(tgse)){
        gvv <- "NA"
        tc <- colnames(tgse)[c]
        for(i in 1:length(gsmval)){
          gsmdati <- as.character(gsub(".*:","",gsmval[i]))
          gsmlabi <- as.character(gsub(":.*","",gsmval[i]))
          # loop over tgse column names
          if(tc %in% gsmlabi){
            gvv <- as.character(gsmdati)
          }
        }
        gvk <- c(gvk,gvv)
      }
      # append new gsm data line
      tgse <- rbind(tgse, gvk)
      message(l)
    }
    tgse.list[[gseid]] <- tgse
  }
  message("gse:",g)
}

#----------
# save data
#----------
# save gse tables list
save(tgse.list, file="geo_gse-atables_list.rda")

# extract as annotation tables
for(t in 1:length(tgse.list)){
  write.csv(tgse.list[[t]], file=paste0(names(tgse.list)[t],"_annotable.csv"))
  message(t)
}

write.csv(gat, file="geo_hmd_22gse_12kgsm.csv")
