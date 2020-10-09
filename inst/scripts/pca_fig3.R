library(data.table)
library(recountmethylation)
library(ggplot2)

# pca data, filenames
pca.all.fn <- "pcadat_allsamples.rda"
pca.noblood.fn <- "pcadat_all-noblood.rda"
pca.nct7.fn <- "pcadat_nct7.rda"
pca.nct6.fn <- "pcadat_nct6-nosperm.rda"

# nct 7 tissues
difid.nct7 <- get(load("gsmid-nct7.rda"))

# get metadata with labels
mdfn <- "metadata_pca-labels.rda"
md <- get(load(mdfn))

# modify labels
md[md$sample.label %in% c("blood", "brain"),]$sample.label <- "other"
md[md$gsm %in% dfid.nct7[dfid.nct7$tissue == "blood",]$gsmid,]$sample.label <- "blood"
md[md$gsm %in% dfid.nct7[dfid.nct7$tissue == "brain",]$gsmid,]$sample.label <- "brain"
table(md$sample.label)


# load feature hashed data, and match to md
fhdat <- fread("fh1000_noobbeta-autosome.table", 
                 header = TRUE, sep = " ", data.table = FALSE)
fhdat <- fhdat[!duplicated(fhdat[,1]),]
rownames(fhdat) <- gsub("\\..*", "", fhdat[,1])
fhdat <- fhdat[,c(2:ncol(fhdat))]
gsm.int <- intersect(md$gsm, rownames(fhdat))
fhdat <- fhdat[rownames(fhdat) %in% gsm.int,]
fhdat <- fhdat[order(match(rownames(fhdat), gsm.int)),]
mdf <- md[md$gsm %in% gsm.int,]
mdf <- mdf[order(match(mdf$gsm, rownames(fhdat))),]
identical(rownames(fhdat), mdf$gsm)

#-----
# pcas
#-----
# all samples
set.seed(2)
pca.all <- prcomp(fhdat)
save(pca.all, file = pca.all.fn)

# all, no blood/leukemia
set.seed(0)
gsm.filt <- md[!md$sample.label %in% c("blood", "leukemia"),]$gsm
fhdat.filt <- fhdat[rownames(fhdat) %in% gsm.filt,]
dim(fhdat.filt)
pca.noblood <- prcomp(fhdat.filt)
plot(pca.noblood$x[,1], pca.noblood$x[,2])
# do rotation so results similar to all samples
plot(pca.noblood$x[,1], -1*pca.noblood$x[,2])
save(pca.noblood, file = pca.noblood.fn)

# nct7 with sperm
set.seed(2)
dfid <- get(load("gsmid-nct7.rda"))
which.samp <- rownames(fhdat) %in% dfid$gsmid
fhdat.filt <- fhdat[which.samp,]
dim(fhdat.filt)
pca.nct7 <- prcomp(fhdat.filt)
save(pca.nct7, file = pca.nct7.fn)

# nct6 no sperm
set.seed(2)
dfid <- dfid[!dfid$tissue == "sperm",]
which.samp <- rownames(fhdat) %in% dfid$gsmid
fhdat.filt <- fhdat[which.samp,]
dim(fhdat.filt)
pca.nct6 <- prcomp(fhdat.filt)
save(pca.nct6, file = pca.nct6.fn)

#---------------------------------
# plots and outliers -- all samples, with blood
#---------------------------------
pca <- get(load(pca.all.fn))
pca.dat <- as.data.frame(pca$x, stringsAsFactors = FALSE)
mdf <- md
identical(mdf$gsm, rownames(pca.dat))

# plots
{
  # get percent var contrib
  sdat <- as.data.frame(summary(pca)$importance, stringsAsFactors = FALSE)
  sperc <- as.numeric(round(100 * sdat[2, c(1:2)], 0))
  cname12 <- paste0(colnames(pca.dat)[1:2], " (", sperc[1:2], "%)")
  
  # add tissue labels
  colvals.facet <- c("blood" = "red", "leukemia" = "purple", "other" = "black")
  collab <- ifelse(mdf$sample.label %in% c("blood", "leukemia"), mdf$sample.label, "other")
  pca.dat$label <- collab
  
  dat1 <- pca.dat[pca.dat$label == "other",]
  dat2 <- pca.dat[pca.dat$label == "blood",]
  dat3 <- pca.dat[pca.dat$label == "leukemia",]
  dat23 <- pca.dat[pca.dat$label %in% c("blood", "leukemia"),]
  colvals <- ifelse(dat23$label == "blood", "red", "purple")
  
  # overplot
  plotname <- "pca-fh1k_overplot_all-bloodlabel.pdf"
  
  pdf(plotname, width = 3, height = 2.8)
  
  ggplot() +
    geom_point(data = dat1, aes(x = PC1, y = PC2), alpha = 0.1) + 
    geom_point(data = dat23, aes(x = PC1, y = PC2), color = colvals) +
    theme_bw() + xlab(cname12[1]) + ylab(cname12[2])
  
  dev.off()
  
  # facet wrap
  plotname <- "pca-fh1k_facet_all-bloodlabel.pdf"
  
  pdf(plotname, width = 4.5, height = 1.8)
  
  ggplot(pca.dat, aes(x = PC1, y = PC2, color = label)) +
    geom_point(alpha = 0.1) + scale_color_manual(values = colvals.facet) +
    theme_bw() + facet_wrap(vars(label)) + xlab(cname12[1]) + ylab(cname12[2]) +
    theme(legend.position = "none")
  
  dev.off()
}

# note outliers
filt1 <- pca.dat[,1] > -10 & mdf$sample.label == "blood"
mdf[filt1,]$gsm
# "GSM2337551", "blood;whole_blood"               
# "GSM2905153", umbilical_cord;blood;cord_blood;stem_cell
# "GSM2465256", blood;whole_blood
# "GSM1954502", lymphatic_system;blood;white_blood_cell;t_cell

table(mdf[filt1,]$gseid)
#   GSE108562;GSE108564   GSE75405;GSE75406   GSE87648;GSE87650            GSE93933 
#         1                   1                   1                   1 

mdf[filt1,]$tissue

filt1 <- pca.dat[,1] > -5 & mdf$sample.label == "blood"
mdf[filt1,]$gsm
# "GSM2337551" "GSM2465256" "GSM3264333" "GSM2465254" "GSM2337461" "GSM2337570"
# whole blood and menstrual blood
table(mdf[filt1,]$gseid)


# studies of menstrual blood (GSE116924, Arai et al 2020),
# inflammatory bowel disease (GSE87648;GSE87650, Ventham et al 2016),
# and cleft lip (GSE93933)

# analysis -- variance
vt1 <- var.test(pca.dat[mdf$sample.label == "blood", 1],
         pca.dat[mdf$sample.label == "leukemia", 1])
vt2 <- var.test(pca.dat[mdf$sample.label == "blood", 2],
                pca.dat[mdf$sample.label == "leukemia", 2])

v1.blood <- var(pca.dat[mdf$sample.label == "blood", 1])
v1.leukemia <- var(pca.dat[mdf$sample.label == "leukemia", 1])
ratio.var1 <- v1.blood/v1.leukemia
ratio.var1

v2.blood <- var(pca.dat[mdf$sample.label == "blood", 2])
v2.leukemia <- var(pca.dat[mdf$sample.label == "leukemia", 2])
ratio.var2 <- v2.blood/v2.leukemia
ratio.var2

# analysis -- means
tt1 <- t.test(pca.dat[mdf$sample.label == "blood", 1],
              pca.dat[mdf$sample.label == "leukemia", 1])
tt2 <- t.test(pca.dat[mdf$sample.label == "blood", 2],
              pca.dat[mdf$sample.label == "leukemia", 2])

#-----------------------------------------
# plots and outliers -- all, without blood, brain labels
#-----------------------------------------
pca <- get(load(pca.noblood.fn))
pca.dat <- as.data.frame(pca$x, stringsAsFactors = FALSE)
# rotate y-axis
pca.dat[,2] <- -1*pca.dat[,2]
# match metadata
mdf <- md[md$gsm %in% rownames(pca.dat),]
mdf <- mdf[order(match(mdf$gsm, rownames(pca.dat))),]
identical(rownames(pca.dat), mdf$gsm)

# plots
{
  # get percent var contrib
  sdat <- as.data.frame(summary(pca)$importance, stringsAsFactors = FALSE)
  sperc <- as.numeric(round(100 * sdat[2, c(1:2)], 0))
  cname12 <- paste0(colnames(pca.dat)[1:2], " (", sperc[1:2], "%)")
  #colnames(pca.dat)[1:2] <- cname12
  
  # add tissue labels
  colvals.facet <- c("brain" = "blue", "brain.tumor" = "cyan4", "other" = "black")
  collab <- ifelse(mdf$sample.label %in% c("brain", "brain.tumor"), mdf$sample.label, "other")
  pca.dat$label <- collab
  
  dat1 <- pca.dat[pca.dat$label == "other",]
  dat2 <- pca.dat[pca.dat$label == "brain",]
  dat3 <- pca.dat[pca.dat$label == "brain.tumor",]
  dat23 <- pca.dat[pca.dat$label %in% c("brain", "brain.tumor"),]
  colvals <- ifelse(dat23$label == "brain", "blue", "cyan4")
  
  # overplot
  plotname <- "pca-fh1k_overplot_all-brain-noblood.pdf"
  
  pdf(plotname, width = 3, height = 2.8)
  
  ggplot() +
    geom_point(data = dat1, aes(x = PC1, y = PC2), alpha = 0.1) + 
    geom_point(data = dat23, aes(x = PC1, y = PC2), color = colvals) +
    theme_bw() + xlab(cname12[1]) + ylab(cname12[2])
  
  dev.off()
  
  # facet wrap
  
  plotname <- "pca-fh1k_facet_all-brain-noblood.pdf"
  
  pdf(plotname, width = 4.5, height = 1.8)
  
  ggplot(pca.dat, aes(x = PC1, y = PC2, color = label)) +
    geom_point(alpha = 0.1) + scale_color_manual(values = colvals.facet) +
    theme_bw() + facet_wrap(vars(label)) +
    xlab(cname12[1]) + ylab(cname12[2]) +
    theme(legend.position = "none")
  
  dev.off()
}

# note outliers
filt1 <- pca.dat[,1] > 0 & pca.dat[,2] < -5 & mdf$sample.label == "brain.tumor"
mdf[filt1,]$gsm 
# "GSM1555032", medulloblastoma,
# "GSM1555033", medulloblastoma, metastatic
# GSM2905381, unknown tumor, metastatic 
# GSM2905379, unknown tumor, metastatic
# GSM2905380, unknown, metastatic
# GSM2905378, unknown, metastatic
# GSM1555030, medulloblastoma

table(mdf[filt1,]$gseid)
# GSE108576 GSE63670;GSE63669 
#     4                 3

# samples from metastatic tumors (6) and medulloblastomas (2)

table(mdf[filt1,]$disease)
# cancer 
#   8
table(mdf[filt1,]$tissue)
# tumor;brain tumor;metastasis;brain 
#   2                      6

# analysis -- variance
vt1 <- var.test(pca.dat[mdf$sample.label == "brain", 1],
                pca.dat[mdf$sample.label == "brain.tumor", 1])
vt2 <- var.test(pca.dat[mdf$sample.label == "brain", 2],
                pca.dat[mdf$sample.label == "brain.tumor", 2])

v1.brain <- var(pca.dat[mdf$sample.label == "brain", 1])
v1.braintumor <- var(pca.dat[mdf$sample.label == "brain.tumor", 1])
ratio.var1 <- v1.brain/v1.braintumor
ratio.var1

v2.brain <- var(pca.dat[mdf$sample.label == "brain", 2])
v2.braintumor <- var(pca.dat[mdf$sample.label == "brain.tumor", 2])
ratio.var2 <- v2.brain/v2.braintumor
ratio.var2

# analysis -- means
tt1 <- t.test(pca.dat[mdf$sample.label == "blood", 1],
              pca.dat[mdf$sample.label == "leukemia", 1])
tt2 <- t.test(pca.dat[mdf$sample.label == "blood", 2],
              pca.dat[mdf$sample.label == "leukemia", 2])

#----------------------------
# plots -- nct7, inc. sperm
#----------------------------
plotname <- "pca-fh1k_nct7.pdf"
# nct7
pca <- get(load(pca.nct7.fn))
pca.dat <- as.data.frame(pca$x, stringsAsFactors = FALSE)
dfid <- get(load("gsmid-nct7.rda"))
dfid <- dfid[order(match(dfid$gsmid, rownames(pca.dat))),]
identical(dfid$gsmid, rownames(pca.dat))

# get percent var contrib
sdat <- as.data.frame(summary(pca)$importance, stringsAsFactors = FALSE)
sperc <- as.numeric(round(100 * sdat[2, c(1:2)], 0))
cname12 <- paste0(colnames(pca.dat)[1:2], " (", sperc[1:2], "%)")

# add col labels
colvals <- c("adipose" = "firebrick", "blood" = "red", 
                   "brain" = "purple", "buccal" = "orange", 
                   "liver" = "forestgreen", "nasal" = "green", "sperm" = "blue")
pca.dat$tissue <- dfid$tissue

# overplot
pdf(plotname, width = 3, height = 3)
ggplot(pca.dat, aes(x = PC1, y = PC2, color = tissue)) +
  geom_point(alpha = 0.4) + theme_bw() + 
  scale_color_manual(values = colvals) +
  xlab(cname12[1]) + ylab(cname12[2]) +
  theme(legend.position = "none")
dev.off()

#----------------------------
# plots -- nct6, without sperm
#----------------------------
plotname <- "pca-fh1k_nct6.pdf"

# nct6
pca <- get(load(pca.nct6.fn))
pca.dat <- as.data.frame(pca$x, stringsAsFactors = FALSE)
dfid <- get(load("gsmid-nct7.rda"))
dfid <- dfid[dfid$gsmid %in% rownames(pca.dat),]
dfid <- dfid[order(match(dfid$gsmid, rownames(pca.dat))),]
identical(dfid$gsmid, rownames(pca.dat))

# get percent var contrib
sdat <- as.data.frame(summary(pca)$importance, stringsAsFactors = FALSE)
sperc <- as.numeric(round(100 * sdat[2, c(1:2)], 0))
cname12 <- paste0(colnames(pca.dat)[1:2], " (", sperc[1:2], "%)")

# add col labels
colvals <- c("adipose" = "firebrick", "blood" = "red", 
             "brain" = "purple", "buccal" = "orange", 
             "liver" = "forestgreen", "nasal" = "green")
pca.dat$tissue <- dfid$tissue

# overplot
pdf(plotname, width = 3, height = 3)
ggplot(pca.dat, aes(x = PC1, y = PC2, color = tissue)) +
  geom_point(alpha = 0.4) + theme_bw() + 
  scale_color_manual(values = colvals) +
  xlab(cname12[1]) + ylab(cname12[2]) +
  theme(legend.position = "none")
dev.off()