#!/usr/bin/env R

# Author: Sean Maden
# Make scatterplot of autosomal DNAm in all tissues without blood or leukemia 
# samples, with labels for brain and brain tumor (fig3b).

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
pca.fn <- "pcadat_all-noblood.rda"
pca <- get(load(file.path(pca.dir, pca.fn)))
pca.dat <- as.data.frame(pca$x, stringsAsFactors = FALSE)
mdf <- mdf[mdf$gsm %in% rownames(pca.dat),]
identical(mdf$gsm, rownames(pca.dat))

#-------------
# prep results
#-------------
# rotate y-axis
pca.dat[,2] <- -1*pca.dat[,2]

# get percent var contrib
sdat <- as.data.frame(summary(pca)$importance, stringsAsFactors = FALSE)
sperc <- as.numeric(round(100 * sdat[2, c(1:2)], 0))
cname12 <- paste0(colnames(pca.dat)[1:2], " (", sperc[1:2], "%)")

# add tissue labels
colvals.facet <- c("brain" = "blue", "brain.tumor" = "cyan4", "other" = "black")
collab <- ifelse(mdf$sample.label %in% c("brain", "brain.tumor"), 
                 mdf$sample.label, "other")
pca.dat$label <- collab

dat1 <- pca.dat[pca.dat$label == "other",]
dat2 <- pca.dat[pca.dat$label == "brain",]
dat3 <- pca.dat[pca.dat$label == "brain.tumor",]
dat23 <- pca.dat[pca.dat$label %in% c("brain", "brain.tumor"),]
colvals <- ifelse(dat23$label == "brain", "blue", "cyan4")

num.samp <- nrow(pca.dat)
plot.title<-paste0("No blood or leukemia samples\n(",num.samp," samples)")
pca.dat$label2 <- pca.dat$label
pca.dat[pca.dat$label == "brain.tumor",]$label2 <- "brain tumor"
collabs <- c("other" = "black", "brain tumor" = "cyan4", "brain" = "blue")

#-----------------
# make scatterplot
#-----------------
fig3b <- ggplot() +
  geom_point(data = dat1, aes(x = PC1, y = PC2), alpha = 0.1) + 
  geom_point(data = dat23, aes(x = PC1, y = PC2), color = colvals) +
  theme_bw() + xlab(cname12[1]) + ylab(cname12[2])

#---------------------
# make manuscript plot
#---------------------
# font sizes
fs.axis.text <- fs.legend.text <- 12; fs.axis.title <- fs.legend.title <- 15

fig3b <- ggplot() +
  geom_point(data = pca.dat[pca.dat$label2 == "other",], 
             aes(x = PC1, y = PC2, color = label2), alpha = 0.2) +
  geom_point(data = pca.dat[pca.dat$label2 %in% c("brain tumor", "brain"),], 
             aes(x = PC1, y = PC2, color = label2)) +
  scale_color_manual(values = collabs) + theme_bw() + 
  xlim(-100, 180) + ylim(-120, 80) +
  ggtitle("All excluding blood and\nleukemia (28,579 samples)") +
  xlab(cname12[1]) + ylab(cname12[2]) + labs(color = "Type") +
  theme(axis.text = element_text(size = fs.axis.text),
        legend.text = element_text(size = fs.legend.text),
        axis.title = element_text(size = fs.axis.title),
        legend.title = element_text(size = fs.legend.title),
        plot.title = element_text(size = fs.axis.title))

# print for manuscript
#pdf("fig3b_pca-fh1k_overplot_all-brain-noblood.pdf", 4.8, 3.5)
#print(fig3b);dev.off()
