#!/usr/bin/env R

# Make scatterplot of autosomal DNAm in all tissues, with labels
# for blood and leukemia samples (fig3a).

library(ggplot2)

#-------------------
# load data, results
#-------------------
pkgname <- "recountmethylationManuscriptSupplement"
pca.dir <- system.file("extdata", "pcadata", package = pkgname) 
# metadata with new labels
md.name <- "metadata_pca-labels.rda"
mdf <- get(load(file.path(pca.dir, "metadata_pca-labels.rda")))
# pca results
pca.fn <- "pcadat_allsamples.rda"
pca <- get(load(file.path(pca.dir, pca.fn)))
pca.dat <- as.data.frame(pca$x, stringsAsFactors = FALSE)
identical(mdf$gsm, rownames(pca.dat))

#-------------
# prep results
#-------------
# get percent var contrib
sdat <- as.data.frame(summary(pca)$importance, stringsAsFactors = FALSE)
sperc <- as.numeric(round(100 * sdat[2, c(1:2)], 0))
cname12 <- paste0(colnames(pca.dat)[1:2], " (", sperc[1:2], "%)")

# add tissue labels
colvals.facet <- c("blood" = "red", "leukemia" = "purple", "other" = "black")
collab <- ifelse(mdf$sample.label %in% c("blood", "leukemia"), mdf$sample.label, "other")
pca.dat$label <- collab
dat1 <- pca.dat[pca.dat$label == "other",]
dat2 <- pca.dat[pca.dat$label == "blood",]
dat3 <- pca.dat[pca.dat$label == "leukemia",]
dat23 <- pca.dat[pca.dat$label %in% c("blood", "leukemia"),]
colvals <- ifelse(dat23$label == "blood", "red", "purple")

#-----------------
# make scatterplot
#-----------------
fig3a <- ggplot() +
  geom_point(data = dat1, aes(x = PC1, y = PC2), alpha = 0.1) + 
  geom_point(data = dat23, aes(x = PC1, y = PC2), color = colvals) +
  theme_bw() + xlab(cname12[1]) + ylab(cname12[2])