#!user/bin/env R

require(jsonlite)

soft2json <- function(gsmsoft_fn_list,gsm_json_destdir){
  for(gsmindex in 1:length(gsmsoft_fn_list)){
    gsmsoft_fn = gsmsoft_fn_list[gsmindex]
    linefile <- read.table(gsmsoft_fn,sep="\n")
    dffile <- as.data.frame(matrix(nrow=1,ncol=nrow(linefile)))
    colnames(dffile) <- gsub(" =.*","",linefile[,1])
    dffile[1,] <- gsub("!.*= ","",linefile[,1])
    jsoni <- toJSON(dffile, pretty=T)
    write(jsoni,file=file.path(gsm_json_destdir,paste0(gsmsoft_fn,".json")))
  }
}

soft2json(gsmsoft_fn_list = commandArgs(T)[1],
          gsm_json_destdir = commandArgs(T)[2])


