#!/usr/bin/env R

# Author: Sean Maden
# Make list object containing study/GSE SOFT-derived metadata. This script 
# reads in filtered JSON-formatted sample/GSM data, derived from GEO GSE 
# SOFT files, and returns as a table in a list. 
# 
# This script takes as input a query object, as generated using the script
# "edirect_query.py", that associates GSM IDs with GSE IDs, which is then 
# used to read in filtered GSM JSON files The input is a filtered table of 
# IDs returned from a call to the script.

require(data.table)
require(rjson)

#---------------------------------
# load data, specify json location
#---------------------------------
# read in a query filt table 
#eq.path <- file.path("recount-methylation-files", "equery", "gsequery_filt.1602784017")
pkgname <- "recountmethylationManuscriptSupplement"
md.dir <- system.file("extdata", "metadata", package = pkgname) 
eq.fn <- "gsequery_filt.1552256509"
eq.path <- file.path(md.dir, eq.fn)

#-----------------------------
# get a list of gsm ids by gse
#-----------------------------
x <- scan(eq.path, what="", sep="\n")
gsel <- list()
for(i in 1:length(x)){
  ssi <- unlist(strsplit(x[[i]]," "))
  gsel[[ssi[1]]] <- ssi[2:length(ssi)]
  message(i)
}

#---------------------------------------------
# get tables list from filtered gsm json files
#---------------------------------------------
dir <- "gsm_json"
dpath <- file.path("recount-methylation-files", dir)
lf.all <- list.files(dpath)
tgse.list <- list()
for(g in 1:length(gsel)){
  lgse <- list(); gseid <- as.character(names(gsel)[g])
  ggse <- gsel[[g]] # gsm list
  # parse json metadata
  for(j in 1:length(ggse)){
    ffilt <- grepl(ggse[j], lf.all) & grepl("\\.filt", lf.all)
    lff <- lf.all[ffilt]
    if(length(lff) > 0){
      fnj <- file.path(dpath, lff)
      lgse[[ggse[j]]] <- fromJSON(paste(readLines(fnj), collapse=""))
    }
  }
  if(length(lgse)>0){
    tcols <- c() # make table columns
    for(l in 1:length(lgse)){
      gsmid <- names(lgse)[l]; gsmval <- unlist(lgse[[l]])
      for(k in 1:length(gsmval)){
        if(grepl(":",gsmval[k])){
          kk <- as.character(gsub(":.*","",gsmval[k]))
          if(!kk=="" & !kk %in% tcols){tcols <- c(tcols, kk)}
        }
      }; message(l)
    }
    tcols <- c("gsm","gse",tcols)
    tgse <- matrix(nrow=0,ncol=length(tcols))
    colnames(tgse) <- tcols
    # coerce to study-specific table
    for(l in 1:length(lgse)){
      gsmid <- names(lgse)[l]; gsmval <- unlist(lgse[[l]]); gvk <- c(gsmid, gseid)
      # loop over gsm values
      for(c in 3:ncol(tgse)){gvv <- "NA"; tc <- colnames(tgse)[c]
        for(i in 1:length(gsmval)){
          gsmdati <- as.character(gsub(".*:","",gsmval[i]))
          gsmlabi <- as.character(gsub(":.*","",gsmval[i]))
          # loop over tgse column names
          if(tc %in% gsmlabi){gvv <- as.character(gsmdati)}
        }; gvk <- c(gvk,gvv)
      }; tgse <- rbind(tgse, gvk); message(l)
    }; tgse.list[[gseid]] <- tgse
  }; message("gse:",g)
}

#----------
# save data
#----------
# save gse tables list
new.fn <- "geo_gse-atables_list"
save(tgse.list, file=paste0(new.fn, ".rda"))

# extract as annotation tables
#for(t in 1:length(tgse.list)){
#  new.fn <- paste0(names(tgse.list)[t],"_annotable.csv")
#  write.csv(tgse.list[[t]], file = new.fn)
#}
