#!/usr/bin/env R

# Make plot of epigenetic vs. mined ages, for panel Fig 1b.

library(recountmethylation)

#----------
# load data
#----------
mdpath <- system.file("extdata", "gsm_metadata",
                      "md_final_hm450k_0-0-1.rda", 
                      package = "recountmethylation")
md <- get(load(mdpath)) # sample metadata

#---------------
# filter samples
#---------------
mdf <- md[!md$age == "valm:NA",]
# get formatted ages
mdf$chron.age <- as.numeric(gsub(";.*", "", gsub("^valm:", "", mdf$age)))
mdf$predage <- as.numeric(mdf$predage)
mdf <- mdf[!is.na(mdf$chron.age),]
mdf <- mdf[!is.na(mdf$predage),]
# make formatted sample type variable
mdf$stype <- as.character(gsub(";.*", "", 
                               gsub("^msraptype:", "", mdf$sampletype)))
mdf <- mdf[!is.na(mdf$stype),]
# make is.cx variable
mdf$is.cx <- ifelse(grepl(".*cancer.*", mdf$disease), TRUE, FALSE)
# get age differences
xdif <- ngsm <- c()
for(g in unique(mdf$gseid)){
  mdff <- mdf[mdf$gseid==g, ]
  xdif <- c(xdif, mean(abs(mdff$chron.age - as.numeric(mdff$predage))))
  ngsm <- c(ngsm, nrow(mdff))}
names(xdif) <- names(ngsm) <- unique(mdf$gseid)
# filter the metadata
filt <- mdf$stype == "tissue" & !mdf$is.cx
filt <- filt & !mdf$gseid %in% names(xdif[xdif > 10])
mdff <- mdf[filt, ]

#----------
# make plot
#----------
fig1b <- ggplot(mdff, aes(x = chron.age, y = predage)) +
  geom_point(size = 1.2, alpha = 0.2) + 
  geom_smooth(method = "lm", size = 1.2) +
  theme_bw() + xlab("Mined chronological age") + 
  ylab("Epigenetic age")

#pdf("fig1b_chron-vs-epigen-age.pdf", 2.7, 2.5)
#print(fig1b); dev.off()