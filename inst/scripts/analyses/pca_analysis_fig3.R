#!/usr/bin/env R

# Author: Sean Maden
# Perform principal component anlayses (PCA) of autosomal DNAm
# across all samples, and within sample subsets. Inputs to this
# script are the large feature-hashed tables of each respective 
# DNAm dataset (access these from 
# https://recount.bio/data/recountmethylation_manuscript_supplement/data/).

library(data.table)
library(recountmethylation)
library(ggplot2)

#----------
# load data
#----------
# results data filenames
pca.all.fn <- "pcadat_allsamples.rda"; pca.noblood.fn <- "pcadat_all-noblood.rda"
pca.nct7.fn <- "pcadat_nct7.rda"; pca.nct6.fn <- "pcadat_nct6-nosperm.rda"

# load feature hashed data, and match to md
fhdat <- fread("fh1000_noobbeta-autosome.table", header = TRUE, sep = " ", data.table = FALSE)
fhdat <- fhdat[!duplicated(fhdat[,1]),]; rownames(fhdat) <- gsub("\\..*", "", fhdat[,1])

# get metadata with labels
pca.dir <- system.file("extdata", "pca", package = pkgname)
mdfn <- "metadata_pca-labels.rda"
md <- get(load(file.path(pca.dir, mdfn)))

# nct 7 tissues
pkgname <- "recountmethylationManuscriptSupplement"
nct7.dir <- system.file("extdata", "nct7", package = pkgname)
fn <- "gsmid-nct7.rda"
difid.nct7 <- get(load(file.path(nct7.dir, fn)))

#----------------
# format metadata
#----------------
# modify labels using nct7 categories
md[md$sample.label %in% c("blood", "brain"),]$sample.label <- "other"
md[md$gsm %in% dfid.nct7[dfid.nct7$tissue == "blood",]$gsmid,]$sample.label <- "blood"
md[md$gsm %in% dfid.nct7[dfid.nct7$tissue == "brain",]$gsmid,]$sample.label <- "brain"
table(md$sample.label)

# match md samples to feature hashed dataset
fhdat <- fhdat[,c(2:ncol(fhdat))]; gsm.int <- intersect(md$gsm, rownames(fhdat))
fhdat <- fhdat[rownames(fhdat) %in% gsm.int,]
fhdat <- fhdat[order(match(rownames(fhdat), gsm.int)),]
mdf <- md[md$gsm %in% gsm.int,]; mdf <- mdf[order(match(mdf$gsm, rownames(fhdat))),]
identical(rownames(fhdat), mdf$gsm)

#------------------------------------
# run pcas, save results data objects
#------------------------------------
set.seed(2)

# all samples
pca.all <- prcomp(fhdat)
pca.all$x <- pca.all$x[,c(1:2)]
save(pca.all, file = pca.all.fn)

# all, no blood/leukemia
gsm.filt <- md[!md$sample.label %in% c("blood", "leukemia"),]$gsm
fhdat.filt <- fhdat[rownames(fhdat) %in% gsm.filt,]; dim(fhdat.filt)
pca.noblood <- prcomp(fhdat.filt)
pca.noblood$x <- pca.noblood$x[,c(1:2)]
save(pca.noblood, file = pca.noblood.fn)

# nct7 with sperm
which.samp <- rownames(fhdat) %in% difid.nct7$gsmid
fhdat.filt <- fhdat[which.samp,]
pca.nct7 <- prcomp(fhdat.filt); save(pca.nct7, file = pca.nct7.fn)

# nct6 no sperm
dfid.nct6 <- difid.nct7[!difid.nct7$tissue == "sperm",]
which.samp <- rownames(fhdat) %in% dfid.nct6$gsmid
fhdat.filt <- fhdat[which.samp,]
pca.nct6 <- prcomp(fhdat.filt); save(pca.nct6, file = pca.nct6.fn)

#---------------
# fig3a analyses
#---------------
# note outliers
pca.dat <- pca.all$x
filt1 <- pca.dat[,1] > -10 & mdf$sample.label == "blood"
mdf[filt1,]$gsm
# "GSM2337551", "blood;whole_blood"               
# "GSM2905153", umbilical_cord;blood;cord_blood;stem_cell
# "GSM2465256", blood;whole_blood
# "GSM1954502", lymphatic_system;blood;white_blood_cell;t_cell

table(mdf[filt1,]$gseid)
#   GSE108562;GSE108564   GSE75405;GSE75406   GSE87648;GSE87650            GSE93933 
#         1                   1                   1                   1 

mdf[filt1,]$tissue

filt1 <- pca.dat[,1] > -5 & mdf$sample.label == "blood"
mdf[filt1,]$gsm
# "GSM2337551" "GSM2465256" "GSM3264333" "GSM2465254" "GSM2337461" "GSM2337570"
# whole blood and menstrual blood

# studies of menstrual blood (GSE116924, Arai et al 2020),
# inflammatory bowel disease (GSE87648;GSE87650, Ventham et al 2016),
# and cleft lip (GSE93933)

# analysis -- variance
vt1 <- var.test(pca.dat[mdf$sample.label == "blood", 1],
                pca.dat[mdf$sample.label == "leukemia", 1])
vt2 <- var.test(pca.dat[mdf$sample.label == "blood", 2],
                pca.dat[mdf$sample.label == "leukemia", 2])

v1.blood <- var(pca.dat[mdf$sample.label == "blood", 1])
v1.leukemia <- var(pca.dat[mdf$sample.label == "leukemia", 1])
ratio.var1 <- v1.blood/v1.leukemia

v2.blood <- var(pca.dat[mdf$sample.label == "blood", 2])
v2.leukemia <- var(pca.dat[mdf$sample.label == "leukemia", 2])
ratio.var2 <- v2.blood/v2.leukemia

#---------------
# fig3b analyses
#---------------
pca <- pca.noblood
pca.dat <- as.data.frame(pca$x, stringsAsFactors = FALSE)

# rotate y-axis
pca.dat[,2] <- -1*pca.dat[,2]

# match metadata
mdf <- md[md$gsm %in% rownames(pca.dat),]
mdf <- mdf[order(match(mdf$gsm, rownames(pca.dat))),]

# note outliers
filt1 <- pca.dat[,1] > 0 & pca.dat[,2] < -5 & 
  mdf$sample.label == "brain.tumor"
mdf[filt1,]$gsm 
# "GSM1555032", medulloblastoma,
# "GSM1555033", medulloblastoma, metastatic
# GSM2905381, unknown tumor, metastatic 
# GSM2905379, unknown tumor, metastatic
# GSM2905380, unknown, metastatic
# GSM2905378, unknown, metastatic
# GSM1555030, medulloblastoma

table(mdf[filt1,]$gseid)
# GSE108576 GSE63670;GSE63669 
#     4                 3

# samples from metastatic tumors (6) and medulloblastomas (2)

table(mdf[filt1,]$disease)
# cancer 
#   8
table(mdf[filt1,]$tissue)
# tumor;brain tumor;metastasis;brain 
#   2                      6

# analysis -- variance
vt1 <- var.test(pca.dat[mdf$sample.label == "brain", 1],
                pca.dat[mdf$sample.label == "brain.tumor", 1])
vt2 <- var.test(pca.dat[mdf$sample.label == "brain", 2],
                pca.dat[mdf$sample.label == "brain.tumor", 2])

v1.brain <- var(pca.dat[mdf$sample.label == "brain", 1])
v1.braintumor <- var(pca.dat[mdf$sample.label == "brain.tumor", 1])
ratio.var1 <- v1.brain/v1.braintumor

v2.brain <- var(pca.dat[mdf$sample.label == "brain", 2])
v2.braintumor <- var(pca.dat[mdf$sample.label == "brain.tumor", 2])
ratio.var2 <- v2.brain/v2.braintumor