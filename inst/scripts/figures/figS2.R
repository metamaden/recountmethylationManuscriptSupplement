#!/usr/bin/env R

# Author: Sean Maden
# Make stacked barplot of MetaSRA-pipeline sample confidences
# and most likely sample type predictions (figS2).

library(ggplot2)

#----------
# load data
#----------
# metadata with new labels
pkgname <- "recountmethylationManuscriptSupplement"
tables.dir <- system.file("extdata", "tables", package = pkgname) 
md.fn <- "table-s1_mdpost_all-gsm-md.csv"
md <- read.csv(file.path(tables.dir, md.fn), header = TRUE)

#--------------
# get plot data
#--------------
mdf <- md[!grepl(".*msraptype:NA.*", md$sampletype),]
stvar <- mdf$sampletype
conf <- as.numeric(gsub(";.*", "", gsub(".*msrapconf:", "", stvar)))
stype <- as.character(gsub(";.*", "", gsub(".*msraptype:", "", stvar)))
bpdf <- matrix(nrow = 0, ncol = 3)
for(i in seq(0, 1, 0.1)){
  for(s in unique(stype)){
    which.gsm <- stype == s & conf >= i
    num.gsm <- nrow(mdf[which.gsm,])
    ms <- matrix(c(s, num.gsm, i), nrow = 1)
    bpdf <- rbind(bpdf, ms)
  }
}

bpdf <- as.data.frame(bpdf, stringsAsFactors = FALSE)
bpdf[,2] <- as.numeric(bpdf[,2])
bpdf[,3] <- as.numeric(bpdf[,3])
colnames(bpdf) <- c("Sample type", "num.gsm", "min.confidence")

#----------------
# get plot object
#----------------
figS2 <- ggplot(bpdf, aes(x = min.confidence, y = num.gsm, fill = `Sample type`)) + 
  geom_bar(stat = "identity") + theme_bw() + xlab("Minimum confidence") +
  ylab("Number of samples") + scale_x_continuous(breaks = seq(0, 1, 0.1))

#pdf("sfig2_bp_stconfcut_msrapst.pdf", 6, 2.2)
#print(figS2);dev.off()
