#!/usr/bin/env R

# Make scatterplot of autosomal DNAm in 7 noncancer tissues
# (fig3c).

library(ggplot2)

#-------------------
# load data, results
#-------------------
pkgname <- "recountmethylationManuscriptSupplement"
pca.dir <- system.file("extdata", "pcadata", package = pkgname) 
nct7.dir <- system.file("extdata", "nct7", package = pkgname) 

# pca results
pca.fn <- "pcadat_nct7.rda"
pca <- get(load(file.path(pca.dir, pca.fn)))
pca.dat <- as.data.frame(pca$x, stringsAsFactors = FALSE)

# nct7
dfid <- get(load(file.path(nct7.dir, "gsmid-nct7.rda")))
dfid <- dfid[order(match(dfid$gsmid, rownames(pca.dat))),]
identical(dfid$gsmid, rownames(pca.dat))

#-------------
# prep results
#-------------
# get percent var contrib
sdat <- as.data.frame(summary(pca)$importance, stringsAsFactors = FALSE)
sperc <- as.numeric(round(100 * sdat[2, c(1:2)], 0))
cname12 <- paste0(colnames(pca.dat)[1:2], " (", sperc[1:2], "%)")

# add col labels
colvals <- c("adipose" = "firebrick", "blood" = "red", 
             "brain" = "purple", "buccal" = "orange", 
             "liver" = "forestgreen", "nasal" = "green", "sperm" = "blue")
pca.dat$tissue <- dfid$tissue

#-----------------
# make scatterplot
#-----------------
fig3c <- ggplot(pca.dat, aes(x = PC1, y = PC2, color = tissue)) +
  geom_point(alpha = 0.4) + theme_bw() + 
  scale_color_manual(values = colvals) +
  xlab(cname12[1]) + ylab(cname12[2]) +
  theme(legend.position = "none")