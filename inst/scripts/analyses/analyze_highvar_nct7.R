#!/usr/bin/env R

# Author: Sean Maden
# Do analyses (binomial tests) of probes with high variance in 7
# noncancer tissues.

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
tname <- "table-s8_14k-highvar_7nct.csv"
ct <- read.csv(file.path(tables.dir, tname))

#---------------
# binomial tests
#---------------
lbm <- list() # results list
padj.max <- 1e-3 # significance cutoff

# genic regions
bm <- matrix(nrow = 0, ncol = 4)
bg.val <- nrow(cga[!cga$UCSC_RefGene_Accession == "",])/nrow(cga)
for(t in unique(ct$tissue)){
  ctt <- ct[ct$tissue == t,]
  ctt.val <- nrow(ctt[!ctt$UCSC_RefGene_Accession == "",])
  btt <- binom.test(ctt.val, n = 2000, p = bg.val)
  mdat <- c(t, btt$p.value, ctt.val/2000, bg.val)
  bm <- rbind(bm, matrix(mdat, nrow = 1))
}
colnames(bm) <- c("tissue", "pval", "fract", "bg")
bm <- as.data.frame(bm, stringsAsFactors = FALSE)
bm$p.adj <- p.adjust(bm[,2], method = "BH")
bm$pos.enrichment <- ifelse(bm$fract > bm$bg, TRUE, FALSE)
bm$sig.pos.enrichment <- ifelse(bm$fract > bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
bm$sig.neg.enrichment <- ifelse(bm$fract < bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
lbm[["genic"]] <- bm

# open sea regions
bm <- matrix(nrow = 0, ncol = 4)
bg.val <- nrow(cga[cga$Relation_to_Island == "OpenSea",])/nrow(cga)
for(t in unique(ct$tissue)){
  ctt <- ct[ct$tissue == t,]
  ctt.val <- nrow(ctt[ctt$Relation_to_Island == "OpenSea",])
  btt <- binom.test(ctt.val, n = 2000, p = bg.val)
  mdat <- c(t, btt$p.value, ctt.val/2000, bg.val)
  bm <- rbind(bm, matrix(mdat, nrow = 1))
}
colnames(bm) <- c("tissue", "pval", "fract", "bg")
bm <- as.data.frame(bm, stringsAsFactors = FALSE)
bm$p.adj <- p.adjust(bm[,2], method = "BH")
bm$pos.enrichment <- ifelse(bm$fract > bm$bg, TRUE, FALSE)
bm$sig.pos.enrichment <- ifelse(bm$fract > bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
bm$sig.neg.enrichment <- ifelse(bm$fract < bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
lbm[["opensea"]] <- bm

# promoter regions
bm <- matrix(nrow = 0, ncol = 4)
bg.val <- nrow(cga[grepl(".*5'.*|.*TSS.*", cga$UCSC_RefGene_Group),])/nrow(cga)
for(t in unique(ct$tissue)){
  ctt <- ct[ct$tissue == t,]
  ctt.val <- nrow(ctt[grepl(".*5'.*|.*TSS.*", ctt$UCSC_RefGene_Group),])
  btt <- binom.test(ctt.val, n = 2000, p = bg.val)
  mdat <- c(t, btt$p.value, ctt.val/2000, bg.val)
  bm <- rbind(bm, matrix(mdat, nrow = 1))
}
colnames(bm) <- c("tissue", "pval", "fract", "bg")
bm <- as.data.frame(bm, stringsAsFactors = FALSE)
bm$p.adj <- p.adjust(bm[,2], method = "BH")
bm$pos.enrichment <- ifelse(bm$fract > bm$bg, TRUE, FALSE)
bm$sig.pos.enrichment <- ifelse(bm$fract > bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
bm$sig.neg.enrichment <- ifelse(bm$fract < bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
lbm[["promoter"]] <- bm

# body regions
bm <- matrix(nrow = 0, ncol = 4)
bg.val <- nrow(cga[!grepl(".*5'.*|.*TSS.*", cga$UCSC_RefGene_Group) &
                     !cga$UCSC_RefGene_Group == "",])/nrow(cga)
for(t in unique(ct$tissue)){
  ctt <- ct[ct$tissue == t,]
  ctt.val <- nrow(ctt[!grepl(".*5'.*|.*TSS.*", ctt$UCSC_RefGene_Group) &
                        !ctt$UCSC_RefGene_Group == "",])
  btt <- binom.test(ctt.val, n = 2000, p = bg.val)
  mdat <- c(t, btt$p.value, ctt.val/2000, bg.val)
  bm <- rbind(bm, matrix(mdat, nrow = 1))
}
colnames(bm) <- c("tissue", "pval", "fract", "bg")
bm <- as.data.frame(bm, stringsAsFactors = FALSE)
bm$p.adj <- p.adjust(bm[,2], method = "BH")
bm$pos.enrichment <- ifelse(bm$fract > bm$bg, TRUE, FALSE)
bm$sig.pos.enrichment <- ifelse(bm$fract > bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
bm$sig.neg.enrichment <- ifelse(bm$fract < bm$bg & 
                                  bm$p.adj < padj.max, TRUE, FALSE)
lbm[["body"]] <- bm

#------------------
# results summaries
#------------------

# pos enrichment for genic regions
lbm$genic[lbm$genic[,7],1]

# neg enrichment for genic regions
lbm$genic[lbm$genic[,8],1] 
# [1] "adipose" "sperm" 

# pos enrichment for opensea regions
lbm$opensea[lbm$opensea[,7], 1] 
# [1] "adipose" "brain"   "buccal"  "liver"   "nasal"   "sperm" 
# neg enrichment for opensea regions
lbm$opensea[lbm$opensea[,8], 1]

# pos enrichment for promoter regions
lbm$promoter[lbm$promoter[,7], 1]
# neg enrichment for promoter regions
lbm$promoter[lbm$promoter[,8], 1]
# [1] "adipose" "brain"   "buccal"  "liver"   "nasal"   "sperm" 

# pos enrichment for body regions
lbm$body[lbm$body[,7], 1]
# [1] "adipose" "brain"   "buccal"  "liver"   "nasal"   "sperm"
# neg enrichment for body regions
lbm$body[lbm$body[,8], 1]
