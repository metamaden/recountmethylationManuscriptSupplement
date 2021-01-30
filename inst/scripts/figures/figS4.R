#!/usr/bin/env R

# Author: Sean Maden
# Get BeadArray control outcomes from quality signals, do PCA and 
# ANOVAs on sample outcomes, then summarize results in component-wise
# stacked barplot screeplot (figS4). 

library(ggplot2)
library(RColorBrewer)
library(recountmethylationManuscriptSupplement)

#----------
# load data
#----------
pkgname <- "recountmethylationManuscriptSupplement"
tables.dir <- system.file("extdata", "tables", package = pkgname)
qcmd.fn <- "table-s2_qcmd-allgsm.csv"
qcmd <- bt <- read.csv(file.path(tables.dir, qcmd.fn), header = TRUE)

bacols <- colnames(qcmd[,grepl("^ba\\..*", colnames(qcmd))])
which.cnames <- c("gsm", "gseid", bacols)
badat <- qcmd[,which.cnames]

#-----------------------
# pca using binary state
#-----------------------
# convert data to numeric
qf <- qcmd[,c(1,3:19)]
colnames(qf) <- gsub("ba\\." , "", colnames(qf))
bt <- bathresh(qf)
btnum <- bt
for(c in 2:ncol(btnum)){btnum[,c] <- ifelse(btnum[,c] == "PASS", 1, 0)}
table(btnum$biotin.stain.red)
colnames(btnum)[2:ncol(btnum)] <- paste0("ba.", colnames(btnum)[2:ncol(btnum)])

# pca
btpc <- prcomp(btnum[,c(2:ncol(btnum))])
btpca.dat <- as.data.frame(btpc$x, stringsAsFactors = FALSE)

# scatterplot
#ggplot(btpca.dat, aes(x = PC1, y = PC2)) +
#  geom_point(alpha = 0.5)

#---------------------
# do variance analysis
#---------------------
vname.dat <- c(gsub("^ba\\.", "", colnames(btnum)[2:ncol(btnum)]), "Residuals")
ssqdat <- matrix(nrow = 0, ncol = length(vname.dat))
colnames(ssqdat) <- vname.dat
bplot <- matrix(nrow = 0, ncol = 3) # barplot df
for(ci in seq(ncol(btpca.dat))){
  pcdat <- btpca.dat[,ci]
  lmpc <- lm(pcdat ~ btnum$ba.restoration.grn + btnum$ba.biotin.stain.red + 
               btnum$ba.biotin.stain.grn + btnum$ba.specificityI.red + btnum$ba.specificityI.grn +
               btnum$ba.specificityII + btnum$ba.extension.red + btnum$ba.extension.grn +
               btnum$ba.hyb.hi.med + btnum$ba.hyb.med.low + btnum$ba.target.removal.1 +
               btnum$ba.target.removal.2 + btnum$ba.bisulfite.conv.I.red +
               btnum$ba.bisulfite.conv.I.grn + btnum$ba.bisulfite.conv.II +
               btnum$ba.nonpolymorphic.red + btnum$ba.nonpolymorphic.grn)
  anpc <- anova(lmpc)
  names.form <- gsub("^btnum\\$ba\\.", "", rownames(anpc))
  nr <- matrix(anpc$`Sum Sq`, nrow = 1)
  names(nr) <- names.form
  nr[!names(nr) %in% vname.dat] <- 0
  # add missing terms
  name.out <- vname.dat[!vname.dat %in% names(nr)]
  nr.append <- rep(0, length(name.out))
  names(nr.append) <- name.out
  nr <- c(nr, nr.append)
  nr <- nr[order(match(names(nr), vname.dat))]
  # barplot
  mdat <- c(names(nr), as.numeric(nr), rep(ci, length(nr)))
  bm <- matrix(mdat, byrow = FALSE, ncol = 3)
  # append new data
  ssqdat <- rbind(ssqdat, nr)
  bplot <- rbind(bplot, bm)
}
rownames(ssqdat) <- paste0("PC", seq(17))
colnames(ssqdat) <- vname.dat
colnames(bplot) <- c("Variable", "ssq", "component")
bplot <- as.data.frame(bplot)
bplot$ssq <- as.numeric(bplot$ssq)

#--------------
# get plot data
#--------------
# get pc var perc
pcvar.sum <- apply(ssqdat, 1, sum)
pcvar.perc <- round(100*pcvar.sum/sum(pcvar.sum), 0)

# modify component labels
ulab <- paste0(gsub("PC", "", names(pcvar.perc)), " (", pcvar.perc, "%)")
bplot$component <- rep(ulab, each = 18)
pcnum <- as.numeric(gsub("PC|\n.*","", bplot$component))
bplot$component <- factor(bplot$component, levels = ulab)

# make stacked barplot
colvect <- c("blue", "red", "green", "purple", "brown", "grey", "pink",
             "firebrick", "cyan", "burlywood", "darkgoldenrod", "darkgreen",
             "darkslategray4", "deeppink", "gray48", "aquamarine", "cadetblue", 
             "chocolate")

#----------
# make plot
#----------
figS4 <- ggplot(bplot, aes(x = component, y = ssq, fill = Variable)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = colvect) + 
  theme_bw() + ylab("Sum of squared variances") + xlab("Component") +
  theme(text = element_text(size=20), 
        axis.text.x = element_text(angle = 90))

#pdf("sfig4_bapca-binnum-thresh.pdf", 8.5, 6)
#print(figS4); dev.off()
