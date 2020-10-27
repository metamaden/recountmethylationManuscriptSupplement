#!/usr/bin/env R

# Get term enrichment by sample BeadArray control performance

library(data.table)
library(recountmethylationManuscriptSupplement)

#----------
# load data
#----------
pkgname <- "recountmethylationManuscriptSupplement"
tables.dir <- system.file("extdata", "tables", package = pkgname)
# load qc data
fn1 <- "table-s2_qcmd-allgsm.csv"
qcmd <- fread(file.path(tables.dir, fn1), sep = ",", data.table = FALSE)
# load metadata
fn2 <- "table-s1_mdpost_all-gsm-md.csv"
md <- fread(file.path(tables.dir, fn2), sep = ",", data.table = FALSE)

#-----------------
# get performances
#-----------------
basignal <- qcmd[,c(1, 3:19)]
colnames(basignal)[2:18] <- gsub("ba\\.", "", colnames(basignal)[2:18])
bat <- bathresh(basignal)

#----------------------------------------------------
# do binomial tests on learned labels, by each metric
#----------------------------------------------------
# get terms where >100 samples below metric threshold
cnv <- colnames(bat)[2:18]; mingsm <- 100; lm <- list()
for(i in 1:length(cnv)){
  mtx <- mdx <- matrix(nrow = 0, ncol = 7)
  c <- cnv[i]; terml <- c(); bai <- bat[bat[,c] == "FAIL", 1]
  if(length(bai) > mingsm){
    terml.tx <- terml.dx <- c(); mdc <- md[md$gsm %in% bai, ]
    terml.tx <- c(terml.tx, unique(unlist(strsplit(mdc$tissue, ";"))))
    terml.dx <- c(terml.dx, unique(unlist(strsplit(mdc$disease, ";"))))
    # na filt
    terml.tx <- terml.tx[!terml.tx == "NA"]; terml.dx <- terml.dx[!terml.dx == "NA"]
    message("working on ", c, " tissue terms...")
    for(t in terml.tx){
      tff <- grepl(paste0(".*", t, ".*"), mdc$tissue)
      bgtf <- grepl(paste0(".*", t, ".*"), md$tissue); bgtt = which(bgtf)
      btest <- binom.test(length(which(tff)), length(tff),
                         p = length(bgtt)/nrow(md))
      # smaple frequencies
      tfreq <- length(which(tff))/nrow(mdc)
      bgfreq <- length(bgtt)/nrow(md)
      mtx <- rbind(mtx, matrix(c(btest$statistic, length(bgtt), 
                                tfreq, bgfreq, ifelse(tfreq > bgfreq, 
                                                      "overenriched", "underenriched"),
                                t, btest$p.value), nrow = 1))
    }
    message("working on ", c, " disease terms...")
    for(d in terml.dx){
      tff <- grepl(paste0(".*", d, ".*"), mdc$disease)
      bgtf <- grepl(paste0(".*", d, ".*"), md$disease); bgtt = which(bgtf)
      btest <- binom.test(length(which(tff)), length(tff),
                         p = length(bgtt)/nrow(md))
      # smaple frequencies
      tfreq <- length(which(tff))/nrow(mdc)
      bgfreq <- length(bgtt)/nrow(md)
      mdx <- rbind(mdx, matrix(c(btest$statistic, length(bgtt), 
                                tfreq, bgfreq, ifelse(tfreq > bgfreq, 
                                                      "overenriched", "underenriched"),
                                d, btest$p.value), nrow = 1))
    }
    colnames(mtx) <- colnames(mdx) <- c("ngsm.set", "ngsm.bg", "ngsm.freq.set", 
                                      "ngsm.freq.bg", "enrichment", "term", "p.unadj")
    lm[[c]] <- list("tissue" = mtx, "disease" = mdx)
  }
  message("finished ", c)
}

#------------------------------------
# get enriched labels, by each metric
#------------------------------------
padj.filt <- 1e-3 # min p.val
enrichment.type <- "overenriched"
loe <- list()
for(c in names(lm)){
  message("working on ", c)
  toe <- doe <- ""; lc <- lm[[c]]; lct <- lc$tissue; lcd <- lc$disease
  lctf <- lct[lct[,5] == enrichment.type, , drop = F]
  lcdf <- lcd[lcd[,5] == enrichment.type, , drop = F]
  if(nrow(lctf) > 1){
    padj <- p.adjust(as.numeric(lctf[, 7]), method = "BH")
    pcond <- padj <= padj.filt; toe <- lctf[pcond, 6]
  }
  if(nrow(lcdf) > 1){
    padj <- p.adjust(as.numeric(lcdf[, 7]), method = "BH")
    pcond <- padj <= padj.filt; doe <- lcdf[pcond, 6]
  }
  loe[[c]] <- list("tissue" = toe, "disease" = doe)
  message("finished with ", c)
}

#--------------------------------
# enrichment tests across metrics
#--------------------------------
min.failed <- 1
ctfailed <- apply(bat[,c(2:18)], 1, function(x){length(x[x=="FAIL"])})
gsm.filt <- bat[ctfailed >= min.failed, 1]
bg.tissue <- unlist(strsplit(md$tissue, ";"))
bg.disease <- unlist(strsplit(md$disease, ";"))
which.test <- which(md$gsm %in% gsm.filt); tot.test <- length(which.test)
test.tissue <- unlist(strsplit(md[which.test,]$tissue, ";"))
test.disease <- unlist(strsplit(md[which.test,]$disease, ";"))
bt <- bd <- matrix(nrow = 0, ncol = 5)
# test tissue terms
for(t in unique(test.tissue)){
  ttfilt <- length(test.tissue[test.tissue == t])
  ttfract <- ttfilt/tot.test
  bgtfilt <- length(bg.tissue[bg.tissue == t])
  bgfract <- bgtfilt/nrow(md)
  btt <- binom.test(ttfilt, tot.test, bgtfilt/nrow(md))
  eterm <- ifelse(ttfract > bgfract, "overenriched", "underenriched")
  mtt <- matrix(c(t, ttfract, bgfract, eterm, btt$p.value), nrow = 1)
  bt <- rbind(bt, mtt)
}
# test disease terms
for(d in unique(test.disease)){
  dtfilt <- length(test.disease[test.disease == d])
  dtfract <- dtfilt/tot.test
  bgdfilt <- length(bg.disease[bg.disease == d])
  bgdfract <- bgdfilt/nrow(md)
  binomtest <- binom.test(dtfilt, tot.test, bgdfract)
  eterm <- ifelse(dtfract > bgdfract, "overenriched", "underenriched")
  mdt <- matrix(c(d, dtfract, bgdfract, eterm, binomtest$p.value), nrow = 1)
  bd <- rbind(bd, mdt)
}
colnames(bt) <- colnames(bd) <- c("label", "test.fract", "bg.fract", "enrichment.type", "test.punadj")

# get significantly overenriched terms
padj.max <- 1e-3; btf <- bt[bt[,4] == "overenriched",]
padj <- p.adjust(as.numeric(btf[,5]), method = "BH")
pcond <- padj <= padj.max
btf[pcond, 1]
# [1] "NA"             "liver"          "prostate"       "mucosa"         "sperm"         
# [6] "adipose"        "visceral"       "adjacent"       "oral"           "nervous_system"
# [11] "umbilical_cord" "cord_blood"     "tongue"         "metastasis" 

padj.max <- 1e-3; bdf <- bd[bd[,4] == "overenriched",]
padj <- p.adjust(as.numeric(bdf[,5]), method = "BH")
pcond <- padj <= padj.max
bdf[pcond, 1]
# [1] "NA"                  "case"                "brain_cancer"        "obese"              
# [5] "prostate_cancer"     "arthritis"           "psoriatic_arthritis"
