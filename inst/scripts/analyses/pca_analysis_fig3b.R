#!/usr/bin/env R

# Analyses of autosomal PCA for brain, brain cancer, and other samples
# (Fig 3b).

library(ggplot2)

#----------
# load data
#----------
pkgname <- "recountmethylationManuscriptSupplement"

# pca results data
pca.dir <- system.file("extdata", "pcadata", package = pkgname)
fn <- "pcadat_all-noblood.rda"
pca.all <- get(load(file.path(pca.dir, fn)))
pca.dat <- pca.all$x
# do rotation on y axis
pca.dat[,2] <- -1*pca.dat[,2]

# metadata
mdfn <- "metadata_pca-labels.rda"
mdf <- get(load(file.path(pca.dir, mdfn)))
mdf <- mdf[mdf$gsm %in% rownames(pca.dat),]
mdf <- mdf[order(match(mdf$gsm, rownames(pca.dat))),]
identical(mdf$gsm, rownames(pca.dat))

#---------------------------
# check distributions -- pc1
#---------------------------
dv <- matrix(nrow = 0, ncol = 2)

dati <- pca.dat[mdf$sample.label == "brain", 1]
dvi <- matrix(c(dati, rep("brain", length(dati))), ncol = 2)
dv <- rbind(dv, dvi)

dati <- pca.dat[mdf$sample.label == "brain.tumor", 1]
dvi <- matrix(c(dati, rep("brain.tumor", length(dati))), ncol = 2)
dv <- rbind(dv, dvi)

dati <- pca.dat[!mdf$sample.label %in% c("brain", "brain.tumor"), 1]
dvi <- matrix(c(dati, rep("other", length(dati))), ncol = 2)
dv <- rbind(dv, dvi)

colnames(dv) <- c("pc1", "samples")
dv <- as.data.frame(dv, stingsAsFactors = FALSE)
dv[,1] <- as.numeric(dv[,1])

# make violin plot
pdf("violin_pc1_groups_fig3b.pdf", 4, 3)
ggplot(dv, aes(x = samples, y = pc1)) + geom_violin()
dev.off()

#---------------------------
# check distributions -- pc2
#---------------------------
dv <- matrix(nrow = 0, ncol = 2)

dati <- pca.dat[mdf$sample.label == "brain", 2]
dvi <- matrix(c(dati, rep("brain", length(dati))), ncol = 2)
dv <- rbind(dv, dvi)

dati <- pca.dat[mdf$sample.label == "brain.tumor", 2]
dvi <- matrix(c(dati, rep("brain.tumor", length(dati))), ncol = 2)
dv <- rbind(dv, dvi)

dati <- pca.dat[!mdf$sample.label %in% c("brain", "brain.tumor"), 2]
dvi <- matrix(c(dati, rep("other", length(dati))), ncol = 2)
dv <- rbind(dv, dvi)

colnames(dv) <- c("pc2", "samples")
dv <- as.data.frame(dv, stingsAsFactors = FALSE)
dv[,1] <- as.numeric(dv[,1])

# make violin plot
pdf("violin_pc2_groups_fig3b.pdf", 4, 3)
ggplot(dv, aes(x = samples, y = pc1)) + geom_violin()
dev.off()

#------------------
# variance analysis
#------------------
# analysis -- variance
vt1 <- var.test(pca.dat[mdf$sample.label == "brain", 1],
                pca.dat[mdf$sample.label == "brain.tumor", 1])
vt2 <- var.test(pca.dat[mdf$sample.label == "brain", 2],
                pca.dat[mdf$sample.label == "brain.tumor", 2])

v1.brain <- var(pca.dat[mdf$sample.label == "brain", 1])
v1.braintumor <- var(pca.dat[mdf$sample.label == "brain.tumor", 1])
ratio.var1 <- v1.braintumor/v1.brain # 16.87976

v2.brain <- var(pca.dat[mdf$sample.label == "brain", 2])
v2.braintumor <- var(pca.dat[mdf$sample.label == "brain.tumor", 2])
ratio.var2 <- v2.braintumor/v2.brain # 25.58693

#---------------
# check outliers
#---------------
filt1 <- pca.dat[,1] > 0 & pca.dat[,2] < -5 & mdf$sample.label == "brain.tumor"

data.frame(mdf[filt1,]$gsm, mdf[filt1,]$gsm_title)
# mdf.filt1....gsm                                       mdf.filt1....gsm_title
# 1       GSM1555032     Primary human brain tumor sample_MDT-AP-2014_methylation
# 2       GSM1555033  Metastasis human brain tumor sample_MDT-AP-2014_methylation
# 3       GSM2905381 Genomic DNA from Uncertain primary tumor brain metastasis-04
# 4       GSM2905379 Genomic DNA from Uncertain primary tumor brain metastasis-02
# 5       GSM2905380 Genomic DNA from Uncertain primary tumor brain metastasis-03
# 6       GSM2905378 Genomic DNA from Uncertain primary tumor brain metastasis-01
# 7       GSM1555030     Primary human brain tumor sample_MDT-AP-2004_methylation
