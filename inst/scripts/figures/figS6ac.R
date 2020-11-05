#!/usr/bin/env R

# Author: Sean Maden
# Make facet plots of data subsets for fig3 PCAs (figS6).

library(ggplot2)

#----------
# load data
#----------
# get path info
pkgname <- "recountmethylationManuscriptSupplement"
pca.dir <- system.file("extdata", "pcadata", package = pkgname) 
tables.dir <- system.file("extdata", "tables", package = pkgname) 
# metadata with new labels
md.name <- "metadata_pca-labels.rda"
md <- get(load(file.path(pca.dir, md.name)))

#--------
# fig S6a
#--------
pca.fn <- "pcadat_allsamples.rda"
pca <- get(load(file.path(pca.dir, pca.fn)))
pca.dat <- as.data.frame(pca$x, stringsAsFactors = FALSE)
# filter metadata
mdf <- md[md$gsm %in% rownames(pca.dat),]
mdf <- mdf[order(match(mdf$gsm, rownames(pca.dat))),]
identical(mdf$gsm, rownames(pca.dat))
# add tissue labels
colvals.facet <- c("blood" = "red", "leukemia" = "purple", "other" = "black")
collab <- ifelse(mdf$sample.label %in% c("blood", "leukemia"), mdf$sample.label, "other")
pca.dat$label <- collab
# get percent var contrib
sdat <- as.data.frame(summary(pca)$importance, stringsAsFactors = FALSE)
sperc <- as.numeric(round(100 * sdat[2, c(1:2)], 0))
cname12 <- paste0(colnames(pca.dat)[1:2], " (", sperc[1:2], "%)")
# make facet wrap plot
figS6a <- ggplot(pca.dat, aes(x = PC1, y = PC2, color = label)) +
  geom_point(alpha = 0.1) + scale_color_manual(values = colvals.facet) +
  theme_bw() + facet_wrap(vars(label)) + xlab(cname12[1]) + ylab(cname12[2]) +
  theme(legend.position = "none")

#--------
# fig S6c
#--------
pca.fn <- "pcadat_all-noblood.rda"
pca <- get(load(file.path(pca.dir, pca.fn)))
pca.dat <- as.data.frame(pca$x, stringsAsFactors = FALSE)
pca.dat$PC2 <- pca.dat$PC2 * -1 # do rotation
# filter metadata
mdf <- md[md$gsm %in% rownames(pca.dat),]
mdf <- mdf[order(match(mdf$gsm, rownames(pca.dat))),]
identical(mdf$gsm, rownames(pca.dat))
# get percent var contrib
sdat <- as.data.frame(summary(pca)$importance, stringsAsFactors = FALSE)
sperc <- as.numeric(round(100 * sdat[2, c(1:2)], 0))
cname12 <- paste0(colnames(pca.dat)[1:2], " (", sperc[1:2], "%)")
#colnames(pca.dat)[1:2] <- cname12
# add tissue labels
colvals.facet <- c("brain" = "blue", "brain.tumor" = "cyan4", "other" = "black")
collab <- ifelse(mdf$sample.label %in% c("brain", "brain.tumor"), mdf$sample.label, "other")
pca.dat$label <- collab
# facet wrap plot
figS6c <- ggplot(pca.dat, aes(x = PC1, y = PC2, color = label)) +
  geom_point(alpha = 0.1) + scale_color_manual(values = colvals.facet) +
  theme_bw() + facet_wrap(vars(label)) +
  xlab(cname12[1]) + ylab(cname12[2]) +
  theme(legend.position = "none")