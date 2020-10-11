#!/usr/bin/env R

# Author: Sean Maden
# Summarize probe and genome region mappings at probes with low variances 
# and mean ranges across 7 noncancer tissues (Table S7).

library(minfiData)

#----------
# load data
#----------
# get minfi-formatted probe annotations
cga <- getAnnotation(get(data("MsetEx")))
# get probe ids
pkgname <- "recountmethylationManuscriptSupplement"
nct7.dir <- system.file("extdata", "nct7", package = pkgname)
id.fn <- "cgids_lowvar-lowmean_nct7.rda"
ptid <- get(load(file.path(nct7.dir, id.fn)))

#---------------
# make new table
#---------------
# get anno and match data
cnames <- c("Name", "chr", "pos", "strand", "Type",	"Islands_Name", "Relation_to_Island",	
            "UCSC_RefGene_Name", "UCSC_RefGene_Accession",	"UCSC_RefGene_Group")
which.cn <- which(colnames(cga) %in% cnames)
cgf <- cga[cga$Name %in% ptid, which.cn]
cgf <- cgf[order(match(rownames(cgf), st$cgid)),]
st <- as.data.frame(cgf, stringsAsFactors = FALSE)

#------
# save
#------
st.fn <- "tableS7_recurr-lowvar-cg-anno_7nct"
save(st, file = paste0(sta.name, ".rda"))
write.csv(st, file = paste0(st.fn, ".csv"), row.names = FALSE)