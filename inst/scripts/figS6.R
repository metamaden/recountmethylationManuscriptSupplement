#!/usr/bin/env R

# PCA results

library(ggplot2)

#--------
# fig S6a
#--------
pca.all.fn <- "pcadat_allsamples.rda"

pca <- get(load(pca.all.fn))
pca.dat <- as.data.frame(pca$x, stringsAsFactors = FALSE)
mdf <- md
identical(mdf$gsm, rownames(pca.dat))

# facet wrap
plotname <- "pca-fh1k_facet_all-bloodlabel.pdf"

pdf(plotname, width = 4.5, height = 1.8)

ggplot(pca.dat, aes(x = PC1, y = PC2, color = label)) +
  geom_point(alpha = 0.1) + scale_color_manual(values = colvals.facet) +
  theme_bw() + facet_wrap(vars(label)) + xlab(cname12[1]) + ylab(cname12[2]) +
  theme(legend.position = "none")

dev.off()

#--------
# fig S6b
#--------
pca <- get(load(pca.noblood.fn))
pca.dat <- as.data.frame(pca$x, stringsAsFactors = FALSE)
# rotate y-axis
pca.dat[,2] <- -1*pca.dat[,2]
# match metadata
mdf <- md[md$gsm %in% rownames(pca.dat),]
mdf <- mdf[order(match(mdf$gsm, rownames(pca.dat))),]
identical(rownames(pca.dat), mdf$gsm)

# facet wrap

plotname <- "pca-fh1k_facet_all-brain-noblood.pdf"

pdf(plotname, width = 4.5, height = 1.8)

ggplot(pca.dat, aes(x = PC1, y = PC2, color = label)) +
  geom_point(alpha = 0.1) + scale_color_manual(values = colvals.facet) +
  theme_bw() + facet_wrap(vars(label)) +
  xlab(cname12[1]) + ylab(cname12[2]) +
  theme(legend.position = "none")

dev.off()