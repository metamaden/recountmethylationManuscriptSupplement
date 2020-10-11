#!/usr/bin/env R

# Author: Sean Maden
# Summarize probes with high tissue-specific autosomal DNAm variance across 7 
# noncancer tissues (Table S8).

library(minfiData)

#----------
# load data
#----------
# get minfi-formatted probe annotations
cga <- getAnnotation(get(data("MsetEx")))
# get probe ids
pkgname <- "recountmethylationManuscriptSupplement"
nct7.dir <- system.file("extdata", "nct7", package = pkgname)
id.fn <- "cgids_highvar_nct7.rda"
ptid <- get(load(file.path(nct7.dir, id.fn)))

#---------------
# make new table
#---------------
# get anno and match data
cnames <- c("Name", "chr", "pos", "strand", "Type",	"Islands_Name", "Relation_to_Island",	
            "UCSC_RefGene_Name", "UCSC_RefGene_Accession",	"UCSC_RefGene_Group")
which.cn <- which(colnames(cga) %in% cnames)
cgf <- cga[cga$Name %in% ptid$cgid, which.cn]
cgf <- cgf[order(match(rownames(cgf), ptid$cgid)),]
st <- cbind(ptid, cgf)

#------
# save
#------
st.fn <- "tableS8_highvar_7nct"
save(st, file = paste0(sta.name, ".rda"))
write.csv(st, file = paste0(st.fn, ".csv"), row.names = FALSE)