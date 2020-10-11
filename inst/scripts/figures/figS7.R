#!/usr/bin/env

# Author: Sean Maden
# Summaries of samples and studies as barplots, 
# and fraction of array probes filtered after

library(ggplot2)

#----------
# load data
#----------
# anova barplot data frame
pkgname <- "recountmethylationManuscriptSupplement"
nct7.dir <- system.file("extdata", "nct7", package = pkgname) 
bpdf.path <- file.path(nct7.dir, "anova-filt_nct.table")
bpdf.anovafilt <- read.table(bpdf.path, sep = " ", header = TRUE)
# noncancer tissues ids
nct7.ids.fn <- "gsmid-nct7.rda"
nct7 <- get(load(file.path(nct7.dir, nct7.ids.fn)))
# samples metadata
tables.dir <- system.file("extdata", "tables", package = pkgname) 
md.fn <- "table-s1_mdpost_all-gsm-md.csv"
md <- read.csv(file.path(tables.dir, md.fn))

#-----------------------------------
# figS7b -- data summaries, barplots
#-----------------------------------
which.samp <- md$gsm %in% nct7$gsmid
mdf <- md[which.samp,]

# count studies
df.gse <- matrix(nrow = 0, ncol = 2)
for(t in unique(nct7$tissue)){
  tgsmv <- nct7[nct7$tissue == t,]$gsmid
  ngse <- length(unique(mdf[mdf$gsm %in% tgsmv,]$gseid))
  df.gse <- rbind(df.gse, matrix(c(t, ngse), nrow = 1))
}
df.gse <- as.data.frame(df.gse, stringsAsFactors = FALSE)
colnames(df.gse) <- c("tissue", "num.gse")
df.gse[,2] <- as.numeric(df.gse[,2])
df.gse$tissue <- factor(df.gse$tissue, 
                        levels = df.gse$tissue[order(df.gse$num.gse)])
df.gse$fillcol <- c("red", "purple", "gold", "forestgreen",
                    "green", "blue", "firebrick")
df.gse <- df.gse[order(df.gse$num.gse),]

# count samples
df.gsm <- matrix(nrow = 0, ncol = 2)
for(t in unique(nct7$tissue)){
  tgsmv <- nct7[nct7$tissue == t,]$gsmid
  ngsm <- length(unique(mdf[mdf$gsm %in% tgsmv,]$gsm))
  df.gsm <- rbind(df.gsm, matrix(c(t, ngsm), nrow = 1))
}
df.gsm <- as.data.frame(df.gsm, stringsAsFactors = FALSE)
colnames(df.gsm) <- c("tissue", "num.gsm")
df.gsm[,2] <- as.numeric(df.gsm[,2])
df.gsm$tissue <- factor(df.gsm$tissue, 
                        levels = df.gsm$tissue[order(df.gsm$num.gsm)])
df.gsm$fillcol <- c("red", "purple", "gold", "forestgreen",
                    "green", "blue", "firebrick")
df.gsm <- df.gsm[order(df.gsm$num.gsm),]

# make plot objects
figS7b.studies <- ggplot(df.gse, aes(x = tissue, y = num.gse, fill = tissue)) + 
  geom_bar(stat = "identity") + theme_bw() + theme(legend.position = "none") +
  scale_fill_manual(values = df.gse$fillcol) + geom_text(aes(label=num.gse), vjust=0)

figS7b.samples <- ggplot(df.gsm, aes(x = tissue, y = num.gsm, fill = tissue)) + 
  geom_bar(stat = "identity") + theme_bw() + theme(legend.position = "none") +
  scale_fill_manual(values = df.gsm$fillcol) + geom_text(aes(label=num.gsm), vjust=0)

#-----------------------------------------
# figS7c -- anova filter, stacked barplots
#-----------------------------------------
# order on removed
bff <- bpdf.anovafilt[bpdf.anovafilt$type == "REMOVED",]
lvl.order <- bff$tissue[order(bff[,1])]
bpdf.anovafilt$tissue <- factor(bpdf.anovafilt$tissue, levels = lvl.order)

# get plot object
figS7c <- ggplot(bpdf.anovafilt, aes(x = tissue, y = nprobes, fill = type)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
