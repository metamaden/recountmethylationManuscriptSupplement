#!/usr/bin/env R

# Metadata term/label enrichment by storage condition quality 
# performance

library(data.table)
library(recountmethylationManuscriptSupplement)

#----------
# load data
#----------
pkgname <- "recountmethylationManuscriptSupplement"
tables.dir <- system.file("extdata", "tables", package = pkgname)
# load qc data
fn1 <- "table-s2_qcmd-allgsm.csv"
qcmd <- fread(file.path(tables.dir, fn1), sep = ",", data.table = FALSE)
# load metadata
fn2 <- "table-s1_mdpost_all-gsm-md.csv"
md <- fread(file.path(tables.dir, fn2), sep = ",", data.table = FALSE)

#---------------------
# get filtered samples
#---------------------

#-------------------------
# label enrichment -- ffpe
#-------------------------

#-----------------------
# label enrichment -- ff
#-----------------------

