#!/usr/bin/env R

# Author: Sean Maden
# Do analyses (binomial tests) of probes with low variance across
# 7 noncancer tissues.

library(minfiData)

#----------
# load data
#----------
# get autosomal dnam annotation
cga <- getAnnotation(get(data("MsetEx")))
cga <- cga[!cga$chr %in% c("chrY", "chrX"),]
# get cgs data
pkgname <- "recountmethylationManuscriptSupplement"
tables.dir <- system.file("extdata", "tables", package = pkgname)
tname <- "table-s7_pantissue-hypovar_7nct.csv"
ct <- read.csv(file.path(tables.dir, tname))

#---------------
# binomial tests
#---------------
lbm <- list() # results list
padj.max <- 1e-3 # significance cutoff
denom <- nrow(ct)

# genic regions
bm <- matrix(nrow = 0, ncol = 3)
bg.val <- nrow(cga[!cga$UCSC_RefGene_Accession == "",])/nrow(cga)
ctt <- ct
ctt.val <- nrow(ctt[!ctt$UCSC_RefGene_Accession == "",])
btt <- binom.test(ctt.val, n = denom, p = bg.val)
mdat <- c(btt$p.value, ctt.val/denom, bg.val)
bm <- rbind(bm, matrix(mdat, nrow = 1))
colnames(bm) <- c("pval", "fract", "bg")
bm <- as.data.frame(bm, stringsAsFactors = FALSE)
bm$p.adj <- p.adjust(bm[,1], method = "BH")
bm$pos.enrichment <- ifelse(bm$fract > bm$bg, TRUE, FALSE)
bm$sig.pos.enrichment <- ifelse(bm$fract > bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
bm$sig.neg.enrichment <- ifelse(bm$fract < bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
lbm[["genic"]] <- bm

# open sea regions
bm <- matrix(nrow = 0, ncol = 3)
bg.val <- nrow(cga[cga$Relation_to_Island == "OpenSea",])/nrow(cga)
ctt <- ct
ctt.val <- nrow(ctt[ctt$Relation_to_Island == "OpenSea",])
btt <- binom.test(ctt.val, n = denom, p = bg.val)
mdat <- c(btt$p.value, ctt.val/denom, bg.val)
bm <- rbind(bm, matrix(mdat, nrow = 1))
colnames(bm) <- c("pval", "fract", "bg")
bm <- as.data.frame(bm, stringsAsFactors = FALSE)
bm$p.adj <- p.adjust(bm[,1], method = "BH")
bm$pos.enrichment <- ifelse(bm$fract > bm$bg, TRUE, FALSE)
bm$sig.pos.enrichment <- ifelse(bm$fract > bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
bm$sig.neg.enrichment <- ifelse(bm$fract < bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
lbm[["opensea"]] <- bm

# promoter regions
bm <- matrix(nrow = 0, ncol = 3)
bg.val <- nrow(cga[grepl(".*5'.*|.*TSS.*", cga$UCSC_RefGene_Group),])/nrow(cga)
ctt <- ct
ctt.val <- nrow(ctt[grepl(".*5'.*|.*TSS.*", ctt$UCSC_RefGene_Group),])
btt <- binom.test(ctt.val, n = denom, p = bg.val)
mdat <- c(btt$p.value, ctt.val/denom, bg.val)
bm <- rbind(bm, matrix(mdat, nrow = 1))
colnames(bm) <- c("pval", "fract", "bg")
bm <- as.data.frame(bm, stringsAsFactors = FALSE)
bm$p.adj <- p.adjust(bm[,1], method = "BH")
bm$pos.enrichment <- ifelse(bm$fract > bm$bg, TRUE, FALSE)
bm$sig.pos.enrichment <- ifelse(bm$fract > bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
bm$sig.neg.enrichment <- ifelse(bm$fract < bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
lbm[["promoter"]] <- bm

# body regions
bm <- matrix(nrow = 0, ncol = 3)
bg.val <- nrow(cga[!grepl(".*5'.*|.*TSS.*", cga$UCSC_RefGene_Group) &
                     !cga$UCSC_RefGene_Group == "",])/nrow(cga)
ctt <- ct
ctt.val <- nrow(ctt[!grepl(".*5'.*|.*TSS.*", ctt$UCSC_RefGene_Group) &
                      !ctt$UCSC_RefGene_Group == "",])
btt <- binom.test(ctt.val, n = denom, p = bg.val)
mdat <- c(btt$p.value, ctt.val/denom, bg.val)
bm <- rbind(bm, matrix(mdat, nrow = 1))
colnames(bm) <- c("pval", "fract", "bg")
bm <- as.data.frame(bm, stringsAsFactors = FALSE)
bm$p.adj <- p.adjust(bm[,1], method = "BH")
bm$pos.enrichment <- ifelse(bm$fract > bm$bg, TRUE, FALSE)
bm$sig.pos.enrichment <- ifelse(bm$fract > bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
bm$sig.neg.enrichment <- ifelse(bm$fract < bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
lbm[["body"]] <- bm

# island (all non-open sea)
bm <- matrix(nrow = 0, ncol = 3)
bg.val <- nrow(cga[!cga$Relation_to_Island == "OpenSea",])/nrow(cga)
ctt <- ct
ctt.val <- nrow(ctt[!ctt$Relation_to_Island == "OpenSea",])
btt <- binom.test(ctt.val, n = denom, p = bg.val)
mdat <- c(btt$p.value, ctt.val/denom, bg.val)
bm <- rbind(bm, matrix(mdat, nrow = 1))
colnames(bm) <- c("pval", "fract", "bg")
bm <- as.data.frame(bm, stringsAsFactors = FALSE)
bm$p.adj <- p.adjust(bm[,1], method = "BH")
bm$pos.enrichment <- ifelse(bm$fract > bm$bg, TRUE, FALSE)
bm$sig.pos.enrichment <- ifelse(bm$fract > bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
bm$sig.neg.enrichment <- ifelse(bm$fract < bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
lbm[["island_all"]] <- bm

# island (restrictive)
bm <- matrix(nrow = 0, ncol = 3)
bg.val <- nrow(cga[cga$Relation_to_Island == "Island",])/nrow(cga)
ctt <- ct
ctt.val <- nrow(ctt[ctt$Relation_to_Island == "Island",])
btt <- binom.test(ctt.val, n = denom, p = bg.val)
mdat <- c(btt$p.value, ctt.val/denom, bg.val)
bm <- rbind(bm, matrix(mdat, nrow = 1))
colnames(bm) <- c("pval", "fract", "bg")
bm <- as.data.frame(bm, stringsAsFactors = FALSE)
bm$p.adj <- p.adjust(bm[,1], method = "BH")
bm$pos.enrichment <- ifelse(bm$fract > bm$bg, TRUE, FALSE)
bm$sig.pos.enrichment <- ifelse(bm$fract > bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
bm$sig.neg.enrichment <- ifelse(bm$fract < bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
lbm[["island"]] <- bm

# island-and-promoter
bm <- matrix(nrow = 0, ncol = 3)
bg.val <- nrow(cga[cga$Relation_to_Island == "Island" & 
                     grepl(".*TSS.*|.*5'.*", cga$UCSC_RefGene_Group),])/nrow(cga)
ctt <- ct
ctt.val <- nrow(ctt[ctt$Relation_to_Island == "Island" &
                      grepl(".*TSS.*|.*5'.*", ctt$UCSC_RefGene_Group),])
btt <- binom.test(ctt.val, n = denom, p = bg.val)
mdat <- c(btt$p.value, ctt.val/denom, bg.val)
bm <- rbind(bm, matrix(mdat, nrow = 1))
colnames(bm) <- c("pval", "fract", "bg")
bm <- as.data.frame(bm, stringsAsFactors = FALSE)
bm$p.adj <- p.adjust(bm[,1], method = "BH")
bm$pos.enrichment <- ifelse(bm$fract > bm$bg, TRUE, FALSE)
bm$sig.pos.enrichment <- ifelse(bm$fract > bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
bm$sig.neg.enrichment <- ifelse(bm$fract < bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
lbm[["island_and_promoter"]] <- bm

# island-and-genic
bm <- matrix(nrow = 0, ncol = 3)
bg.val <- nrow(cga[cga$Relation_to_Island == "Island" & 
                     !cga$UCSC_RefGene_Group == "",])/nrow(cga)
ctt <- ct
ctt.val <- nrow(ctt[ctt$Relation_to_Island == "Island" &
                      !cga$UCSC_RefGene_Group == "",])
btt <- binom.test(ctt.val, n = denom, p = bg.val)
mdat <- c(btt$p.value, ctt.val/denom, bg.val)
bm <- rbind(bm, matrix(mdat, nrow = 1))
colnames(bm) <- c("pval", "fract", "bg")
bm <- as.data.frame(bm, stringsAsFactors = FALSE)
bm$p.adj <- p.adjust(bm[,1], method = "BH")
bm$pos.enrichment <- ifelse(bm$fract > bm$bg, TRUE, FALSE)
bm$sig.pos.enrichment <- ifelse(bm$fract > bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
bm$sig.neg.enrichment <- ifelse(bm$fract < bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
lbm[["island_and_genic"]] <- bm
