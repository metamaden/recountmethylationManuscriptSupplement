#!/usr/bin/env R

# Author: Sean Maden
# Make violin plots of high variance probes in 7 noncancer tissues 
# (fig4a and fig4b).

library(ggplot2)

#----------
# load data
#----------
pkgname <- "recountmethylationManuscriptSupplement"
nct7.dir <- system.file("extdata", "nct7", package = pkgname) 
tx.fn <- "ltxfilt-dnamsummary_highvar_nct7.rda"
txfiltl <- get(load(file.path(nct7.dir, tx.fn)))
txnames <- names(txfiltl)

#-----------------------
# prep the plot datasets
#-----------------------
bpdf.mean <- matrix(nrow = 0, ncol = 3)
bpdf.var <- matrix(nrow = 0, ncol = 3)

# tx colors
tl <- c("blood", "buccal", "brain", "sperm", "nasal", 
        "adipose", "liver" ) # index
tc <- c("red", "orange", "purple", "blue", "green", 
        "brown", "forestgreen") # color

# get summaries by tissue
for(t in 1:length(txnames)){
  tx <- txnames[t]; datt <- txfiltl[[tx]]
  bpdf.mean <- rbind(bpdf.mean, data.frame(datt$mean, tc[tl == tx], 
                                           tx, stringsAsFactors = F))
  bpdf.var <- rbind(bpdf.var, data.frame(datt$var, tc[tl == tx], 
                                         tx, stringsAsFactors = F))
}

# format plot data
bpdf.mean <- as.data.frame(bpdf.mean, stringsAsFactors = F)
bpdf.var <- as.data.frame(bpdf.var, stringsAsFactors = F)
bpdf.mean[,1] <- as.numeric(as.character(bpdf.mean[,1]))
bpdf.var[,1] <- as.numeric(as.character(bpdf.var[,1]))
colnames(bpdf.mean) <- c("mean", "col", "tissue")
colnames(bpdf.var) <- c("var", "col", "tissue")
bpdf.mean$tissue <- factor(bpdf.mean$tissue, levels = unique(bpdf.mean$tissue))
ordert <- order(match(tl, levels(bpdf.mean$tissue)))
bpdf.var$tissue <- factor(bpdf.var$tissue, levels = unique(bpdf.var$tissue))
ordert <- order(match(tl, levels(bpdf.var$tissue)))

#-------------------
# make plot objects
#-------------------
# means violin plots
fig4a <- ggplot(bpdf.mean, aes(x = tissue, y = mean, fill = tissue)) + 
  geom_violin(trim = F, show.legend = F) +
  scale_fill_manual(values = tc[ordert]) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("Tissue") + ylab("Beta-value mean") +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, 0.5))

# variances violin plots
fig4b <- ggplot(bpdf.var, aes(x = tissue, y = var, fill = tissue)) + 
  geom_violin(trim = F, show.legend = F) +
  scale_fill_manual(values = tc[ordert]) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  xlab("Tissue") + ylab("Beta-value variance") +
  scale_y_continuous(limits = c(0, 0.12), breaks = seq(0, 0.12, 0.06))

# print for manuscript
#pdf("fig4a_vp-mean_hivar_7nct.pdf", 3.2, 2.3)
#print(fig4a); dev.off()
#pdf("fig4b_vp-var_hivar_7nct.pdf", 3.2, 2.3)
#print(fig4b); dev.off()