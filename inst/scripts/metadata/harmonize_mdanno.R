# GSM sample metadata (rsheet) table inc. GSE and MetaSRA-pipeline mapped terms
load("~/scratch/analysis_scratch/analysis_final_files/rm_samplemetadata_mainfiles/current/rsheet_36kGSM.rda")
# GSE-wise semi-manual annotation of sample metadata from GEO soft files
load("~/scratch/analysis_scratch/analysis_final_files/rm_samplemetadata_mainfiles/current/geoanno_478gse.rda")
# MetaSRA-pipeline common mapped terms definitions
load("~/scratch/analysis_scratch/analysis_final_files/rm_samplemetadata_mainfiles/current/md-toi-df_for-val.rda") 

# data prep
dfg <- as.data.frame(gat, stringsAsFactors = F)
dfg$gsmid <- dfg$gsm

# check repeated sample data
tg <- as.data.frame(table(dfg$gsm))
tg <- tg[rev(order(tg[,2])),]
head(tg)

dfi <- dfg[dfg$gsm=="GSM3440790",]
identical(as.character(dfi[1,]), 
          as.character(dfi[2,]))

tgf <- tg[tg[,2]>1,]
sv <- c()
for(t in 1:nrow(tgf)){
  dfi <- dfg[dfg$gsm==tgf[t,1],]; dfi <- dfi[,c(1,3:ncol(dfi))]
  sv <- c(sv,identical(as.character(dfi[1,]), 
            as.character(dfi[2,])))
  message(t)
}
table(sv)
sv
#FALSE  TRUE 
#4397  6312
dfg[dfg$gsm==tgf[which(sv==FALSE)[1],1],]

# notes: 
# data apparently equivalent across duplicates
# use naive duplicates filter

dfgf <- dfg[!duplicated(dfg$gsm),]
dim(dfgf)
# [1] 33362     9

#--------------------------------------
# GSM sample titles from GEO soft data
#--------------------------------------
# iterate on gsm_json_filt files
# scrape just the sample title from filtered soft files
require(rjson)
sfp <- "gsm_json_filt"
lf <- list.files(sfp)
lff <- lf[grepl(".*\\.json\\.filt$",lf)]

gsmtitledf <- data.frame(matrix(c("",""),nrow=0,ncol=2))
for(j in 1:length(lff)){
  gsmtitledf <- rbind(gsmtitledf,
                      data.frame(matrix(c(gsub("\\..*","",lff[j]),
                                          as.character(unlist(fromJSON(paste(readLines(paste0(sfp, "/", lff[j])),collapse="")))["!Sample_title"])),
                                        nrow=1,ncol=2)))
  message(j)
}

tdf <- as.data.frame(gsmtitledf, stringsAsFactors=F)
tdf$gsmid <- as.character(tdf[,1])
save(tdf, file="gsm_jsontitledf.rda")

#-----------------------
# harmonize gsm metadata
#-----------------------
# 0.1 grab all GSM IDs and for mdf
gsm.v <- unique(c(tdf$gsmid,
         rs$gsmid,
         dfgf$gsm))
mdf <- data.frame(gsm=gsm.v, stringsAsFactors = F)
dim(mdf) # 37481 1

# 0.2 prepare all metadata tables to append
dfll <- list(tdf, rs, dfgf)
dfnl <- list()
for(d in 1:length(dfll)){
  odf <- dfll[[d]]
  ngsm <- mdf$gsm[!mdf$gsm %in% odf$gsmid] # outersect gsm ids
  # new df
  ndf <- as.data.frame(matrix(rep("NA",length(ngsm)*ncol(odf)),
                nrow=length(ngsm)))
  colnames(ndf) <- colnames(odf)
  ndf[,which(colnames(odf)=="gsmid")] <- ngsm
  oodf <- rbind(odf, ndf)
  oodf <- oodf[order(match(oodf$gsmid, mdf$gsm)),] # match order with main table
  dfnl[[d]] <- oodf
  message(d,";",identical(oodf$gsmid, mdf$gsm))
}

# GSE manual annotation data
colnames(dfnl[[3]]); 
# [1] "gsm"               "gse"               "sample_type"       "disease_state"     "gender"           
# [6] "age"               "anatomic_location" "misc"              "gsmid"  

# idat/array file metadata from rsheet
colnames(dfnl[[2]])
# [1] "gsmid"            "gseid"            "idats_fn"         "msrapmd_fn"       "msrapmd_flatjson"
# [6] "SENTRIX_ID"       "ARRAY_ID"         "Basename"

# 1. Add annotations 
mdf$gseid <- dfnl[[2]]$gseid
mdf$gsm_title <- as.character(dfnl[[1]]$X2) # sample title
mdf$sample_type <- dfnl[[3]]$sample_type
mdf$disease_state <- dfnl[[3]]$disease_state
mdf$gender <- dfnl[[3]]$gender
mdf$age <- dfnl[[3]]$age 
mdf$anatomic_location <- dfnl[[3]]$anatomic_location
mdf$misc <- dfnl[[3]]$misc
mdf$idats_fn <- dfnl[[2]]$idats_fn
mdf$msrapmd_fn <- dfnl[[2]]$msrapmd_fn
mdf$sentrix_id <- dfnl[[2]]$SENTRIX_ID
mdf$array_id <- dfnl[[2]]$ARRAY_ID
mdf$basename <- dfnl[[2]]$Basename

save(mdf, file="gsm_geomd_hdf.rda")
write.csv(mdf, file="gsm_geomd_harmonizedtable.csv")

#---------------------------
# append processed metadata
#---------------------------
# 1. Append binarized msrap terms for most frequent terms

# 2.1 Make gross sample labels for pca




# 2.2 Identify largest studies lacking gross labels, for further manual annotation

#----------------------------------------------
# Manual anno. largest GSEs lacking gross labels
#----------------------------------------------




#-----------------------------
# append new methyl. metadata
#-----------------------------
# 1. Sex prediction

# 2. Horvath imputed biological age

# 3. Idat batch metadata


