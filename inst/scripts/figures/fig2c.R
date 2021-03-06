#!/usr/bin/env R

# Shows how to make barplot of BeadArray control metric failure
# percentages by storage condition type (Fig 2C).

library(ggplot2)
library(data.table)
library(recountmethylationManuscriptSupplement)

#----------
# load data
#----------
# tables dir
pkgname <- "recountmethylationManuscriptSupplement"
tables.dir <- system.file("extdata", "tables", package = pkgname)
# qc signals
fn <- "table-s2_qcmd-allgsm.csv"
qcmd <- fread(file.path(tables.dir, fn), sep = ",", data.table = FALSE)
# metadata table
tname <- "table-s1_mdpost_all-gsm-md.csv"
md <- read.csv(file.path(tables.dir, tname), header = TRUE)

#---------------------------
# get storage condition data
#---------------------------
mf <- md[!is.na(md$storage),];gsmv <- mf$gsm
gsmv.frozen <- mf[grepl("frozen", mf$storage),]$gsm
gsmv.ffpe <- mf[grepl("FFPE", mf$storage),]$gsm
qf <- qcmd[qcmd$gsm %in% gsmv,];qf <- qf[order(match(qf$gsm, gsmv)),]
if(!identical(qf$gsm, gsmv)){stop("error matching qf to gsmv!")}
rownames(qf) <- qf$gsm;bf <- qf[,grepl("^ba\\..*", colnames(qf))]
colnames(bf) <- gsub("^ba\\.", "", colnames(bf))

#----------------------------------------
# get the control performances, plot data
#----------------------------------------
tf <- bathresh(bf);bpdf <- matrix(nrow = 0, ncol = 3)
lab.frozen <- "Fresh frozen";lab.ffpe <- "FFPE"
for(c in colnames(tf)){
  cd <- tf[rownames(tf) %in% gsmv.frozen, c]
  mfail.frozen <- matrix(c(length(cd[cd=="FAIL"])/length(gsmv.frozen), 
                           c, lab.frozen), nrow = 1)
  cd <- tf[rownames(tf) %in% gsmv.ffpe, c]
  mfail.ffpe <- matrix(c(length(cd[cd=="FAIL"])/length(gsmv.ffpe), 
                         c, lab.ffpe), nrow = 1)
  if(!(mfail.ffpe[1] == 0 & mfail.frozen[1] == 0)){
    bpdf <- rbind(bpdf, rbind(mfail.frozen, mfail.ffpe))}}
colnames(bpdf) <- c("perc.fail", "ba.control", "storage.cond")
bpdf <- as.data.frame(bpdf, stringsAsFactors = FALSE)
bpdf[,1] <- as.numeric(bpdf[,1]);bpdf[,1] <- 100 * bpdf[,1]
# order on total failed samples
uvar <- unique(as.character(bpdf[,2]))
sumv <- unlist(lapply(uvar, function(x){sum(bpdf[bpdf[,2] == x, 1])}))
# varval <- bpdf[bpdf$storage.cond == lab.ffpe,1]
bpdf[,2] <- factor(bpdf[,2], levels = uvar[rev(order(sumv))])
bpdf[,3] <- factor(bpdf[,3], levels = c(lab.frozen, lab.ffpe))

#----------
# make plot
#----------
fig2c <- ggplot(bpdf, aes(x = ba.control, y = perc.fail, 
                          fill = storage.cond)) + 
  geom_bar(stat = "identity", position = "dodge") + theme_bw() +
  scale_fill_manual(values = c("purple", "orange")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  xlab("BeadArray control") + ylab("Failed\nsamples (%)") +
  labs(fill = "Storage\ncondition")

# print for manuscript
#pdf("fig2c_bp-gsmfail_ffpe-ff.pdf", 5, 3)
#print(fig2c); dev.off()
  
