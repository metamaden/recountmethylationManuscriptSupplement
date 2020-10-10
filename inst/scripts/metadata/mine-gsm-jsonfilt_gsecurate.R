# mine gsm json filt data for gse curation
library(rjson)
gset <- read.csv("large_gse_info.csv")
dn <- "gsm_jsonfilt_gsecurate"

gselist <- gset$Var1
tgse.list <- list()

for(g in 1:length(gselist)){
  lgse <- list()
  gseid <- as.character(gset$Var1[g])
  ggse <- gmap[gmap$gseid==gset$Var1[g],]$gsmid
  # parse json metadata
  for(j in 1:length(ggse)){
    lgse[[ggse[j]]] <- fromJSON(paste(readLines(paste0(dn,
                                                       "/",ggse[j],
                                                       ".json.filt")), 
                                      collapse="")
    )
    message(j)
  }
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
  message("gse:",g)
}

save(tgse.list, file="gsecurate_mdtables_list.rda")

