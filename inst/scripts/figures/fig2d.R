#!/usr/bin/env R

# Author: Sean Maden
# Make panel Fig 2C, study-level quality assessment outcomes

library(ComplexHeatmap);library(circlize)
library(recountmethylationManuscriptSupplement)

#----------
# load data
#----------
pkgname <- "recountmethylationManuscriptSupplement"
tables.dir <- system.file("extdata", "tables", package = pkgname)
fn <- "table-s5_qcfreq_gsewise.csv"; 
st.agg <- read.csv(file.path(tables.dir, fn))
# tables.dir <- ""; st.agg <- read.csv(fn)

#----------
# make plot
#----------
# format colnames
which.cv <- grep("meth|unmeth", colnames(st.agg))
cvf <- colnames(st.agg)[which.cv]
colnames(st.agg)[which.cv] <- gsub("_", " ", gsub("\\.", "<", cvf))

# get variables, data matrices
gsesize <- st.agg$ngsm; names(gsesize) <- st.agg$gseid
hm.baf <- st.agg[,3:19]; hm.qcf <- st.agg[,20:25]

# reorder on size
identical(rownames(hm.baf), rownames(hm.qcf)); ro <- rev(order(gsesize))

# get color palettes
hmdat <- hm.baf; breaks <- seq(min(hmdat), max(hmdat), 0.01)
hmcol1 <- colorRamp2(breaks, 
                     colorRampPalette(c("blue",
                                        "orange",
                                        "red"))(n=length(breaks)))
hmdat <- hm.qcf; breaks <- seq(min(hmdat), max(hmdat), 0.01)
hmcol2 <- colorRamp2(breaks, 
                     colorRampPalette(c("black",
                                        "forestgreen",
                                        "green"))(n=length(breaks)))

# get aggregate assessment
qc.stf.gseanno <- ifelse(st.agg$qc.stf.filt == "above", TRUE, FALSE)
ba.stf.gseanno <- ifelse(st.agg$ba.stf.filt == "above", TRUE, FALSE)
agg.stf <- ifelse(ba.stf.gseanno | qc.stf.gseanno, "above", "below")

#----------------------------------------------------
# make composite hm object, studies with fst > thresh
#----------------------------------------------------
# filter studies
gsefilt <- agg.stf == "above"
hm.baf.subset <- hm.baf[gsefilt,];hm.qcf.subset <- hm.qcf[gsefilt,]
gsesize.subset <- gsesize[gsefilt];ro <- rev(order(gsesize.subset))

# get image objects
h1 <- Heatmap(hm.baf.subset[ro,], name = "BeadArray\nS.T. Freq.",
             col = hmcol1, cluster_rows = F, cluster_columns = T,
             show_row_dend = F, show_column_dend = F,
             show_row_names = F)

h2 <- Heatmap(hm.qcf.subset[ro,], name = "Signal\nS.T. Freq.",
             col = hmcol2,
             cluster_rows = F, cluster_columns = T,
             show_row_dend = F, show_column_dend = F,
             show_row_names = F)

hm.ranno <- rowAnnotation(`Log2 Study Size\n(Num. GSM)` = 
                            anno_barplot(log2(gsesize.subset[ro]), 
                                         width = unit(4, "cm"),
                                         xlab = "GSE Size",
                                         col = "red"))

# make composite plot object
hmlist <- h1 + h2 + hm.ranno

#-----------------------------------------
# filter studies >10 samples, above thresh
#-----------------------------------------
# filter studies
min.gsm <- 11; gsefilt <- agg.stf == "above" & st.agg$ngsm >= min.gsm
hm.baf.subset <- hm.baf[gsefilt,];hm.qcf.subset <- hm.qcf[gsefilt,]
gsesize.subset <- gsesize[gsefilt];ro <- rev(order(gsesize.subset))

# get image objects
h1.legend.title <- list(title=expression(italic("BeadArray f"["st"])))
h2.legend.title <- list(title=expression(italic("Signal f"["st"])))
h1 <- Heatmap(hm.baf.subset[ro,], heatmap_legend_param = h1.legend.title,
              col = hmcol1, cluster_rows = F, cluster_columns = T,
              show_row_dend = F, show_column_dend = F, show_row_names = F)

h2 <- Heatmap(hm.qcf.subset[ro,], heatmap_legend_param = h2.legend.title,
              col = hmcol2, cluster_rows = F, cluster_columns = T, 
              show_row_dend = F, show_column_dend = F, show_row_names = F)

hm.ranno <- rowAnnotation("Log2 study size\n(number of samples)" = 
                            anno_barplot(log2(gsesize.subset[ro]), 
                                         width = unit(4, "cm"),
                                         xlab = "GSE Size", col = "red"))
names(hm.ranno) <- expression("log"[2])

# make composite plot object
hmlist <- h1 + h2 + hm.ranno

#---------------------
# save manuscript plot
#---------------------
row.title <- paste0(paste(rep(" ", 28),collapse = ""), 
                    "Study (", length(gsesize.subset)," studies)")

# print for manuscript
#pdf("hmcomp_subset-thresh-min-11gsm_gse-stfreq.pdf", 6.8, 3.5)
#draw(hmlist, row_title = row.title, heatmap_legend_side = "right", 
#     annotation_legend_side = "right"); dev.off()

# draw for vignette
draw(hmlist, row_title = row.title, heatmap_legend_side = "right", 
     annotation_legend_side = "right")