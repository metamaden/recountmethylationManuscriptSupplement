#!/usr/bin/env R

# Make study-wise quality controls table S5. Calculates
# fi for 17 BeadArray controls, methylated and unmethylated signals,
# then applies criteria i (fi, failing at least 1 BA control, > 60%) 
# and criteria ii (fi, below 11 for both signals, > 60%).

library(recountmethylationManuscriptSupplement)

# name of new supp table
supp.table.name <- "table-s5_qcfreq_gsewise.csv"

#----------
# load data
#----------
pkgname <- "recountmethylationManuscriptSupplement"
tables.dir <- system.file("extdata", "tables", package = pkgname)
fn <- "table-s2_qcmd-allgsm.csv"
qcmd <- read.csv(file.path(tables.dir, fn))
bam = qcmd[,grepl("^ba\\..*|^gsm$", colnames(qcmd))]
colnames(bam) = gsub("ba\\.", "", colnames(bam))

#------------------
# get study metrics
#------------------
# study sizes
gsev = unique(qcmd$gseid)
gsesize = c()
for(g in gsev){gsesize = c(gsesize, length(which(qcmd$gseid==g)))}
st.agg = data.frame(gseid = gsev, ngsm = gsesize) # agg for supp table

# BA metrics, sub-threshold freq
gm <- matrix(nrow = 0, ncol = 17)
bat <- bathresh(bam)
rownames(bat) <- bat$gsm
bat <- bat[,c(2:ncol(bat))]
for(g in gsev){
  ngsm <- nrow(qcmd[qcmd$gseid==g,])
  which.gsm <- rownames(bat) %in% qcmd[qcmd$gseid==g,]$gsm
  bag <- apply(bat[which.gsm,], 2, function(x){
    length(x[x=="FAIL"])/ngsm
  })
  gm <- rbind(gm, matrix(c(bag), nrow = 1))
  # message(g)
}
colnames(gm) <- colnames(bat)
rownames(gm) <- gsev
class(gm) <- "numeric"
st.agg <- cbind(st.agg, hm.baf) # agg for supp table

# signal sub 10, 11, 12
gsev <- unique(qcmd$gseid)
bat <- qcmd[,c("meth.l2med", "unmeth.l2med")]
rownames(bat) <- qcmd$gsm
lvlv <- c(10, 11, 12)
gm <- matrix(nrow = 0, ncol = 2*length(lvlv))
for(g in gsev){
  ngsm <- nrow(qcmd[qcmd$gseid==g,])
  which.gsm <- rownames(bat) %in% qcmd[qcmd$gseid==g,]$gsm
  bag <- c()
  for(l in lvlv){
    bag <- c(bag, apply(bat[which.gsm,], 2, function(x){
        length(x[x < l])/ngsm
    }))
  }
  gm <- rbind(gm, matrix(c(bag), nrow = 1))
  #message(g)
}
colnames(gm) <- c("meth_<_10", "unmeth_<_10",
                 "meth_<_11", "unmeth_<_11",
                 "meth_<_12", "unmeth_<_12")
rownames(gm) <- gsev
class(gm) <- "numeric"
st.agg <- cbind(st.agg, hm.qcf) # agg for supp table
st.agg$gseid <- rownames(st.agg)

#-------------------------
# add assessments i and ii
#-------------------------
min.fi <- 0.6

# criterium i, beadarray metrics
hm.baf <- st.agg[,3:19]
thresh = 0.6
ba.stf.gseanno = c()
for(g in 1:nrow(hm.baf)){
  ba.stf.gseanno = c(ba.stf.gseanno,
                     ifelse(length(which(hm.baf[g,] > thresh)) > 0, 
                            "above", "below"))
}

# criterium ii, meth and unmeth signal
hm.qcf <- st.agg[,20:25]
thresh = 0.6
qc.stf.gseanno = c()
which.col = c(3, 4)
for(g in 1:nrow(hm.qcf)){
  qc.stf.gseanno = c(qc.stf.gseanno,
                     ifelse(length(which(hm.qcf[g, which.col] > thresh)) == 2, 
                            "above", "below"))
}

st.agg$mu.stf.filt <- qc.stf.gseanno
st.agg$ba.stf.filt <- ba.stf.gseanno

#---------------------
# save the supp table
#---------------------
write.csv(st.agg, file = supp.table.name)