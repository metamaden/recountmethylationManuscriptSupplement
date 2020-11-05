#!/usr/bin/env R

# Do differential methylome variation analyses of cancers and normals

library(recountmethylationManuscriptSupplement)
library(HDF5Array)
library(ggplot2)

#----------
# load data
#----------
# get set ids
pkgname <- "recountmethylationManuscriptSupplement"
pca.dir <- system.file("extdata", "pcadata", package = pkgname)
md.fn <- "metadata_pca-labels.rda"
md <- get(load(file.path(pca.dir, md.fn)))
# load data
db.path <- "remethdb-h5se_gr_0-0-1_1590090412"
grset <- loadHDF5SummarizedExperiment(db.path)

#-----------------------
# blood versus leukemias
#-----------------------
# get filtered datasets
labv <- c("blood", "leukemia"); mdf <-  md[md$sample.label %in% labv,]
gsmv <- mdf$gsm[mdf$gsm %in% grset$gsm]; mdf <- mdf[mdf$gsm %in% gsmv,]
grf <- grset[,grset$gsm %in% gsmv]; grf <- grf[,order(match(grf$gsm, mdf$gsm))]
identical(grf$gsm, mdf$gsm)
pData(grf) <- DataFrame(mdf)
# get autosomal dnam
anno <- getAnnotation(grf); annof <- anno[!anno$chr %in% c("chrY", "chrX"),]
grf <- grf[rownames(grf) %in% rownames(annof),]
# do var analysis
mv <- matrix(nrow = 0, ncol = 2)
bv <- getblocks(nrow(grf), bsize = 100000)
#bi <- getBeta(grf)
for(b in bv){
  betai <- as.matrix(getBeta(grf[unlist(b),]))
  mi <- matrix(nrow = nrow(betai), ncol = 2)
  for(r in seq(nrow(betai))){
    x <- betai[r,]
    v1 <- var(x[grf$sample.label == "blood"])
    v2 <- var(x[grf$sample.label == "leukemia"])
    mi[r,] <- c(v1, v2)
    message(r)
  }
  mv <- rbind(mv, mi)
}
rownames(mv) <- rownames(grf)
colnames(mv) <- c("var.blood", "var.leukemia")
# assess results
nrow(mv[mv[,1] < mv[,2],]) # 311127
nrow(mv[mv[,1] < mv[,2],])/nrow(mv) # 0.6565745
save(mv, file = "m-var_blood-vs-leukemia.rda")

#-------------------------
# brain versus brain tumor
#-------------------------
# get filtered datasets
labv <- c("brain", "brain.tumor"); mdf <-  md[md$sample.label %in% labv,]
gsmv <- mdf$gsm[mdf$gsm %in% grset$gsm]; mdf <- mdf[mdf$gsm %in% gsmv,]
grf <- grset[,grset$gsm %in% gsmv]; grf <- grf[,order(match(grf$gsm, mdf$gsm))]
identical(grf$gsm, mdf$gsm)
pData(grf) <- DataFrame(mdf)
# get autosomal dnam
anno <- getAnnotation(grf); annof <- anno[!anno$chr %in% c("chrY", "chrX"),]
grf <- grf[rownames(grf) %in% rownames(annof),]
# do var analysis
mv <- matrix(nrow = 0, ncol = 2)
bv <- getblocks(nrow(grf), bsize = 100000)
for(b in bv){
  betai <- as.matrix(getBeta(grf[unlist(b),]))
  mi <- matrix(nrow = nrow(betai), ncol = 2)
  for(r in seq(nrow(betai))){
    x <- betai[r,]
    v1 <- var(x[grf$sample.label == "brain"])
    v2 <- var(x[grf$sample.label == "brain.tumor"])
    mi[r,] <- c(v1, v2)
    message(r)
  }
  mv <- rbind(mv, mi)
}
rownames(mv) <- rownames(grf)
colnames(mv) <- c("var.brain", "var.braintumor")
# assess results
nrow(mv[mv[,1] < mv[,2],]) # 444304
nrow(mv[mv[,1] < mv[,2],])/nrow(mv) # 0.9376192
save(mv, file = "m-var_brain-vs-braintumor.rda")