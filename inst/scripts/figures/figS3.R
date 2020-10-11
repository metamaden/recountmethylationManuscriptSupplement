#!/usr/bin/env R

# Author: Sean Maden
# Make signal plots of BeadArray controls.

#----------
# load data
#----------
pkgname <- "recountmethylationManuscriptSupplement"
tables.dir <- system.file("extdata", "tables", package = pkgname)
qcmd.fn <- "table-s2_qcmd-allgsm.csv"
qcmd <- bt <- read.csv(file.path(tables.dir, qcmd.fn), header = TRUE)
qcmd <- qcmd[,c(3:19)]; colnames(qcmd) <- gsub("ba\\.", "", colnames(qcmd))

# get the control thresholds
dft <- data.frame(restoration.grn = 0, biotin.stain.red = 5, biotin.stain.grn = 5,
                 specificityI.red = 1, specificityI.grn = 1, specificityII = 1,
                 extension.red = 5, extension.grn = 5, hyb.hi.med = 1, 
                 hyb.med.low = 1, target.removal.1 = 1, target.removal.2 = 1, 
                 bisulfite.conv.I.red = 1, bisulfite.conv.I.grn = 1, 
                 bisulfite.conv.II = 1, nonpolymorphic.red = 5, nonpolymorphic.grn = 5, 
                 stringsAsFactors = F)

# specify controls of interest
which.bim <- c("biotin.stain.red", "biotin.stain.grn", "bisulfite.conv.I.red", 
               "bisulfite.conv.I.grn", "bisulfite.conv.II")

# get signals for controls of interest
qcf <- qcmd[,colnames(qcmd) %in% which.bim]

#--------------------
# make composite plot
#--------------------
par(mfcol = c(2, 3), 
    mar = c(2.5, 2.5, 2.5, 2), 
    oma = c(3, 3, 1, 1))
for(c in colnames(qcf)){
  datc <- qcf[, c]
  if(c %in% c("restoration.grn", "biotin.stain.red", "biotin.stain.grn")){
    datc <- datc[datc < quantile(datc, seq(0, 1, 0.1))[8]]
    plot(density(datc), main = paste0(c, "\n(trimmed)"), xlab = "", ylab = "")
  } else{
    plot(density(datc), main = c, xlab = "", ylab = "")
  }
  abline(v = as.numeric(dft[,c]), lwd = 2, col = "blue")
}
mtext(paste0("Signal"), side = 1, outer = T)
mtext("Density", side = 2, outer = T)
dev.off()