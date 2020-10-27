#!/usr/bin/env R

# Analyses of autosomal PCA for blood, leukemia, and other samples
# (Fig 3a).

library(ggplot2)

#----------
# load data
#----------
pkgname <- "recountmethylationManuscriptSupplement"

# pca results data
pca.dir <- system.file("extdata", "pcadata", package = pkgname)
fn <- "pcadat_allsamples.rda"
pca.all <- get(load(file.path(pca.dir, fn)))
pca.dat <- pca.all$x

# metadata
mdfn <- "metadata_pca-labels.rda"
mdf <- get(load(file.path(pca.dir, mdfn)))
mdf <- mdf[order(match(mdf$gsm, rownames(pca.dat))),]
identical(mdf$gsm, rownames(pca.dat))

#---------------------------
# check distributions -- pc1
#---------------------------
dv <- matrix(nrow = 0, ncol = 2)

dati <- pca.dat[mdf$sample.label == "leukemia", 1]
dvi <- matrix(c(dati, rep("leukemia", length(dati))), ncol = 2)
dv <- rbind(dv, dvi)

dati <- pca.dat[mdf$sample.label == "blood", 1]
dvi <- matrix(c(dati, rep("blood", length(dati))), ncol = 2)
dv <- rbind(dv, dvi)

dati <- pca.dat[!mdf$sample.label %in% c("leukemia", "blood"), 1]
dvi <- matrix(c(dati, rep("other", length(dati))), ncol = 2)
dv <- rbind(dv, dvi)

colnames(dv) <- c("pc1", "samples")
dv <- as.data.frame(dv, stingsAsFactors = FALSE)
dv[,1] <- as.numeric(dv[,1])

# violin plot
pdf("violin_pc1_groups_fig3a.pdf", 4, 3)
ggplot(dv, aes(x = samples, y = pc1)) + geom_violin()
dev.off()

#---------------------------
# check distributions -- pc2
#---------------------------
dv <- matrix(nrow = 0, ncol = 2)

dati <- pca.dat[mdf$sample.label == "leukemia", 2]
dvi <- matrix(c(dati, rep("leukemia", length(dati))), ncol = 2)
dv <- rbind(dv, dvi)

dati <- pca.dat[mdf$sample.label == "blood", 2]
dvi <- matrix(c(dati, rep("blood", length(dati))), ncol = 2)
dv <- rbind(dv, dvi)

dati <- pca.dat[!mdf$sample.label %in% c("leukemia", "blood"), 2]
dvi <- matrix(c(dati, rep("other", length(dati))), ncol = 2)
dv <- rbind(dv, dvi)

colnames(dv) <- c("pc2", "samples")
dv <- as.data.frame(dv, stingsAsFactors = FALSE)
dv[,1] <- as.numeric(dv[,1])

# violin plot
pdf("violin_pc2_groups_fig3a.pdf", 4, 3)
ggplot(dv, aes(x = samples, y = pc1)) + geom_violin()
dev.off()

#------------------
# variance analyses
#------------------
v1.blood <- var(pca.dat[mdf$sample.label == "blood", 1])
v1.leukemia <- var(pca.dat[mdf$sample.label == "leukemia", 1])
ratio.var1 <- v1.leukemia/v1.blood # 1.245173

v2.blood <- var(pca.dat[mdf$sample.label == "blood", 2])
v2.leukemia <- var(pca.dat[mdf$sample.label == "leukemia", 2])
ratio.var2 <- v2.leukemia/v2.blood # 6.260585

# ftests for differential variances
vt1 <- var.test(pca.dat[mdf$sample.label == "blood", 1],
                pca.dat[mdf$sample.label == "leukemia", 1])
vt2 <- var.test(pca.dat[mdf$sample.label == "blood", 2],
                pca.dat[mdf$sample.label == "leukemia", 2])

#---------------
# check outliers
#---------------
filt1 <- pca.dat[,1] > -10 & mdf$sample.label == "blood"
data.frame(mdf[filt1,]$gsm, mdf[filt1,]$gsm_title)

# mdf.filt1....gsm                            mdf.filt1....gsm_title
# 1       GSM2337551                              2095CD [Whole blood]
# 2       GSM2465256 genomic DNA from whole-blood DNA from sample  C3b
# 3       GSM3264333                                 Primary-Edom22_P7
# 4       GSM2465254 genomic DNA from whole-blood DNA from sample  C1b
# 5       GSM2337461                               0068H [Whole blood]
# 6       GSM2337570                              0057HC [Whole blood]
