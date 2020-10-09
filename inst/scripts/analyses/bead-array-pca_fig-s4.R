#!/usr/bin/env R

# Author: Sean Maden

# Performs PCA on binary BeadArray control outcomes,
# and plots screeplot ANOVA outcomes.

library(ggplot2)
library(RColorBrewer)

fn.binary <- "qcmd_allmetrics_allgsm.rda"
bn.signal <- "qcmd_allmetrics_allgsm.rda"
bt <- get(load(fn.binary))
qcmd <- get(load(bn.signal))

bacols <- colnames(qcmd[,grepl("^ba\\..*", colnames(qcmd))])
which.cnames <- c("gsm", "gseid", bacols)
badat <- qcmd[,which.cnames]

# font size for plots
theme_set(
  theme_classic(base_size = 25)
)

#---------------------
# pca using raw signal
#---------------------
bapca <- prcomp(badat[,c(3:ncol(badat))])
bapca.dat <- as.data.frame(bapca$x, stringsAsFactors = FALSE)

# scatterplot
ggplot(bapca.dat, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.5)

# variance analysis
ssqdat <- matrix(nrow = 0, ncol = 17)
bplot <- matrix(nrow = 0, ncol = 3) # barplot df
for(ci in seq(ncol(bapca.dat))){
  pcdat <- bapca.dat[,ci]
  lmpc <- lm(pcdat ~ badat$ba.restoration.grn + badat$ba.biotin.stain.red + 
               badat$ba.biotin.stain.grn + badat$ba.specificityI.red + badat$ba.specificityI.grn +
               badat$ba.specificityII + badat$ba.extension.red + badat$ba.extension.grn +
               badat$ba.hyb.hi.med + badat$ba.hyb.med.low + badat$ba.target.removal.1 +
               badat$ba.target.removal.2 + badat$ba.bisulfite.conv.I.red +
               badat$ba.bisulfite.conv.I.grn + badat$ba.bisulfite.conv.II +
               badat$ba.nonpolymorphic.red + badat$ba.nonpolymorphic.grn)
  anpc <- anova(lmpc)
  ssqdat <- rbind(ssqdat, matrix(anpc$`Sum Sq`[1:17], nrow = 1))
  # barplot
  bm <- matrix(c(rownames(anpc)[1:17], 
                 anpc$`Sum Sq`[1:17],
                 rep(paste0("PC", ci), 17)), 
               byrow = FALSE, ncol = 3)
  bplot <- rbind(bplot, bm)
}
rownames(ssqdat) <- paste0("PC", seq(17))
colnames(ssqdat) <- rownames(anpc)[1:17]
colnames(bplot) <- c("metric", "ssq", "pc")
bplot <- as.data.frame(bplot)
bplot$ssq <- as.numeric(bplot$ssq)

# stacked barplot (x-axis component, y-axis variance)
bplot$component <- as.numeric(gsub("PC", "", bplot$pc))
bplot$metric <- gsub("^badat\\$ba\\." , "", bplot$metric)

jpeg("bapca-signal.jpg", 8, 5, units = "in", res = 400)
ggplot(bplot, aes(x = component, y = ssq, fill = metric)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1))
dev.off()


bpf <- bplot[grepl("PC1$|PC2$", bplot$pc),]
ggplot(bpf, aes(x = pc, y = ssq, fill = metric)) +
  geom_bar(stat = "identity", position = "stack")

# metric contribution by component
ssqdf <- as.data.frame(ssqdat, stringsAsFactors = FALSE)
colnames(ssqdf) <- gsub("^badat\\$", "", colnames(ssqdf))
ssqdf$ba.biotin.stain.red[1]/sum(ssqdf[1,]) # 0.9961851
ssqdf$ba.biotin.stain.grn[2]/sum(ssqdf[2,]) # 0.9965472
# component percent of total var
sum(ssqdf[1,])/sum(apply(ssqdf, 1, sum)) # 0.8538052
sum(ssqdf[2,])/sum(apply(ssqdf, 1, sum)) # 0.1461813

# stacked barplot (color by component)
ssq.compsum <-apply(ssqdf, 1, sum)
bpdf.comp <- data.frame(var = as.numeric(ssq.compsum),
                        component = names(ssq.compsum),
                        lab = rep("component", length(ssq.compsum)),
                        stringsAsFactors = FALSE)
comp.order <- order(as.numeric(gsub("PC", "", bpdf.comp$component)))
bpdf.comp$component <- factor(bpdf.comp$component,
                              levels = bpdf.comp$component[comp.order])
ggplot(bpdf.comp, aes(x = lab, y = var, fill = component)) + 
  geom_bar(stat = "identity", position = "stack")

#-----------------------
# pca using binary state
#-----------------------
# convert data to numeric
btnum <- bt
for(c in 2:ncol(btnum)){btnum[,c] <- ifelse(btnum[,c] == "PASS", 1, 0)}
table(btnum$biotin.stain.red)
colnames(btnum)[2:ncol(btnum)] <- paste0("ba.", colnames(btnum)[2:ncol(btnum)])

# pca
btpc <- prcomp(btnum[,c(2:ncol(btnum))])
btpca.dat <- as.data.frame(btpc$x, stringsAsFactors = FALSE)

# scatterplot
ggplot(btpca.dat, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.5)

# var analysis
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
jpeg("bapca-binnum-thresh.jpg", 9, 6, units = "in", res = 400)
ggplot(bplot, aes(x = component, y = ssq, fill = Variable)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = colvect) + 
  theme_bw() + ylab("Sum of Squared Variances") + xlab("Component") +
  theme(text = element_text(size=20), 
        axis.text.x = element_text(angle = 90))
dev.off()

lc <- list()
pcv <- c()
for(i in unique(bplot$component)){
  bf <- bplot[bplot[,3] == i,]
  bf[,2] <- bf[,2]/sum(bf[,2])
  bf <- bf[bf[,1] %in% c("biotin.stain.grn", "biotin.stain.red", 
                         "bisulfite.conv.I.red", "nonpolymorphic.grn",
                         "nonpolymorphic.red"),]
  message(i, ": ", sum(bf[,2]))
  pcv <- c(pcv, sum(bf[,2]))
  lc[[i]] <- paste0(bf[,1], ": ", 100*bf[,2])
}





