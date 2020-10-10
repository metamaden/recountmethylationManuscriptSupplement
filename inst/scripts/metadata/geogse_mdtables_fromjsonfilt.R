# build geo md table by study id, or for all studies

eqf <- "gsequery_filt.1552256509"
# Read in the data
x <- scan(eqf, what="", sep="\n")
gsel <- list()
for(i in 1:length(x)){
  ssi <- unlist(strsplit(x[[i]]," "))
  gsel[[ssi[1]]] <- ssi[2:length(ssi)]
  message(i)
}

dn <- "gsm_json_filt"
gsel <- list()

#----------------
# specify gse id
#----------------
#gse.target <- "GSE49032;GSE49031"
#gse.target <- "GSE49032"
gse.target <- "GSE49031"
lgse <- list()
#gseid <- as.character(names(gsel)[g])
gseid.check <- which(names(gsel)==gse.target)

is.null(gseid.check) # check whether id in names 

# ggse <- gsel[[g]]
ggse <- gsel[[gseid.check]]

length(ggse) # check >0 gsm ids

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
  gsel[[gseid]] <- tgse
}

# define function
get.mdtable <- function(gseid="GSE49031", dn="gsm_json_filt", parse.eq=FALSE, return.eq=FALSE, eqfn=""){
  # Form a GSE-specific metadata table (mdtable) from soft-extracted GSM json files
  # Arguments:
  # gseid (str, required) : valid GSE, must be present in provided equery file
  # dn (str, required): Directory name for GSM filtered JSON metadata files, from GEO soft file
  # parse.eq (bool, required): Whether to parse an edirect query (eq) file
  # return.eq (bool, required): Whether to return parsed eqfile
  # eqfn (str, optional) : Edirect query filename
  # Returns:
  # List or Metadata table: List (parsed eq file and mdtable), or table (prepared mdtable), 
  #   see argument 'return.eq'

  # option to parse gse query file
  if(parse.eq & !eqfn==""){
    message("Parsing equery file...")
    # eqf <- "gsequery_filt.1552256509"
    # Read in the data
    x <- scan(eqfn, what="", sep="\n")
    gsel <- list()
    for(i in 1:length(x)){
      ssi <- unlist(strsplit(x[[i]]," "))
      gsel[[ssi[1]]] <- ssi[2:length(ssi)]
      message(i)
    }
    message("Finished parsing equery file. Continuing...")
  } else{
    message("Skipping equery file parse step...")
  }
  
  lgse <- list()
  gseid.index <- which(names(gsel)==gseid)
  # check if provided gse id in parsed eq file
  if(is.null(gseid.index)){
    message("GSE ID not identified in equery file. Returning...")
    if(return.eq){
      message("Returning parsed equery file...")
      return(gsel)
    } else{
      return()
    }
  }
  
  ggse <- gsel[[gseid.index]]
  # check >0 gsm ids
  if(length(ggse)>0){
    message("Sample GSM IDs identified. Continuing...")
  } else{
    message("No sample GSM IDs found. Returning...")
    if(return.eq){
      message("Returning parsed equery file...")
      return(gsel)
    } else{
      return()
    }
  } 
  
  # parse json filt metadata
  require(rjson)
  for(j in 1:length(ggse)){
    fnj <- paste0("./", dn, "/", ggse[j], ".json.filt")
    if(file.exists(fnj)){
      lgse[[ggse[j]]] <- fromJSON(paste(readLines(fnj), collapse=""))
    }
    message(j)
  }
  if(length(lgse)>0){
    # make table columns
    message("making mdtable columns...")
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
        message("finished gsm val ",k)
      }
      message("finished gsm (index) ",l)
    }
    tcols <- c("gsm","gse",tcols)
    tgse <- matrix(nrow=0,ncol=length(tcols))
    colnames(tgse) <- tcols
    # coerce to study-specific table
    message("corecing data to study mdtable...")
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
      message("finished gsm (index) ", l)
    }
  }
  
  message("finished mdtable prep. Returning...")
  if(return.eq){
    message("Returning GSE table with parsed equery file...")
    lr <- list(eq.parse=gsel, gseid=tgse)
    return(lr)
  } else{
    message("Returning  GSE mdtable...")
    return(tgse)
  }
}
save(get.mdtable,file="get_mdtable_function.rda")

xt <- get.mdtable(gseid="GSE49031")

#--------------------
# all gse, md tables
#--------------------
# form the gse tables list
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