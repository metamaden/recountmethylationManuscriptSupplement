#!/usr/bin/env R

# aggregate metadata

library(data.table)

#-----------------
# cell predictions
#-----------------
# Note: summarize and aggregate 3 reps
cp1 = as.data.frame(fread("cellpred_table", sep = ' ', header = T))
cp2 = as.data.frame(fread("cellpred_table_rep2", sep = ' ', header = T))
cp3 = as.data.frame(fread("cellpred_table_rep3", sep = ' ', header = T))

rownames(cp1) <- cp1[,1]; cp1 <- cp1[,c(2:ncol(cp1))]
rownames(cp2) <- cp2[,1]; cp2 <- cp2[,c(2:ncol(cp2))]
rownames(cp3) <- cp3[,1]; cp3 <- cp3[,c(2:ncol(cp3))]

cp2 <- cp2[order(match(rownames(cp2), rownames(cp1))),]
cp3 <- cp3[order(match(rownames(cp3), rownames(cp1))),]

identical(rownames(cp1), rownames(cp2))
identical(rownames(cp1), rownames(cp3))

par(mfrow=c(3,2))
for(c in 1:6){
  plot(cp1[,c], cp2[,c], main=colnames(cp1)[c])
}


nt = matrix(nrow=0,ncol=6)
for(r in 1:nrow(cp1)){
  gcc <- c()
  for(c in 1:6){
    gcc <- c(gcc, var(c(cp1[r,c], cp2[r,c], cp3[r,c])))
  }
  nt <- rbind(nt, gcc)
  message(r)
}
colnames(nt) <- colnames(cp1)[1:6]
rownames(nt) <- rownames(cp1)
colMeans(nt)
# CD8T         CD4T           NK        Bcell         Mono         Gran 
# 3.238845e-05 1.181664e-05 9.754414e-06 6.010645e-06 9.463399e-06 1.633350e-05 

apply(nt, 2, median)
# CD8T         CD4T           NK        Bcell         Mono         Gran 
# 1.170579e-05 4.123941e-06 3.615454e-06 2.061127e-06 3.439118e-06 5.498687e-06

#-----------------
# sex predictions
#-----------------

xy = as.data.frame(fread("xypred_table", sep = ' ', header = T))

hadf <- data.frame(id = rownames(ha),
                   horvath.predage = ha[,1],
                   stringsAsFactors = F)

# get gsm ids
hadf$gsm <- gsub(".*hlink\\.|_.*", "", hadf$id)
xy$gsm <- gsub(".*hlink\\.|_.*", "", xy$id)
cp$gsm <- gsub(".*hlink\\.|_.*", "", cp$id)

# collate new data
xy <- xy[order(match(xy$gsm, hadf$gsm)),]
identical(xy$gsm, hadf$gsm)
gsmna <- xy$gsm[!xy$gsm %in% cp$gsm]
cp <- rbind(cp, data.frame(id = rep("NA", length(gsmna)),
                           CD8T = rep("NA", length(gsmna)),
                           CD4T = rep("NA", length(gsmna)),
                           NK = rep("NA", length(gsmna)),
                           Bcell = rep("NA", length(gsmna)),
                           Mono = rep("NA", length(gsmna)),
                           Gran = rep("NA", length(gsmna)),
                           gsm = gsmna, stringsAsFactors = F))

cp <- cp[order(match(cp$gsm, xy$gsm)),]
identical(cp$gsm, xy$gsm)
ndf <- data.frame(gsm = xy$gsm,
                  predage = hadf$horvath.predage,
                  predsex = xy$predsex,
                  predcell.CD8T = cp$CD8T,
                  predcell.CD4T = cp$CD4T,
                  predcell.NK = cp$NK,
                  predcell.Bcell = cp$Bcell,
                  predcell.Mono = cp$Mono,
                  predcell.Gran = cp$Gran,
                  stringsAsFactors = F)

# match and append new data
ld = ldff
ndf <- ndf[order(match(ndf$gsm, ld$gsm)),]
# append missing samples
gsmna <- ld$gsm[!ld$gsm %in% ndf$gsm]
ndf <- rbind(ndf, data.frame(gsm = gsmna,
                             predage = rep("NA", length(gsmna)),
                             predsex = rep("NA", length(gsmna)),
                             predcell.CD8T = rep("NA", length(gsmna)),
                             predcell.CD4T = rep("NA", length(gsmna)),
                             predcell.NK = rep("NA", length(gsmna)),
                             predcell.Bcell = rep("NA", length(gsmna)),
                             predcell.Mono = rep("NA", length(gsmna)),
                             predcell.Gran = rep("NA", length(gsmna)),
                             stringsAsFactors = F))
ndf <- ndf[order(match(ndf$gsm, ld$gsm)),]
identical(ndf$gsm, ld$gsm)
rownames(ld) <- ld$gsm
rownames(ndf) <- ndf$gsm

ld <- cbind(ld, ndf)

table(ld$gender, ld$predsex)

mdf <- ld
save(mdf, file="mdfinal_anno_pred.rda")

