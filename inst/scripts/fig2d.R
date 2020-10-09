#!/usr/bin/env R

# Make panel Fig 2C, study-level quality assessment outcomes

library(ComplexHeatmap)
library(circlize)
library(recountmethylationManuscriptSupplement)

#----------
# load data
#----------
pkgname <- "recountmethylationManuscriptSupplement"
tables.dir <- system.file("extdata", "tables", package = pkgname)
fn <- "table-s5_qcfreq_gsewise.csv"
st.agg <- read.csv(file.path(tables.dir, fn))

#----------
# make plot
#----------
gsesize <- st.agg$ngsm
names(gsesize) <- st.agg$gseid
hm.baf <- st.agg[,3:19]
hm.qcf <- st.agg[,20:25]

# reorder on size
identical(rownames(hm.baf), rownames(hm.qcf))
ro <- rev(order(gsesize))

# get color palettes
hmdat <- hm.baf
breaks <- seq(min(hmdat), max(hmdat), 0.01)
hmcol1 <- colorRamp2(breaks, 
                     colorRampPalette(c("blue","orange","red"))(n=length(breaks)))
hmdat <- hm.qcf
breaks <- seq(min(hmdat), max(hmdat), 0.01)
hmcol2 <- colorRamp2(breaks, 
                     colorRampPalette(c("black","green","white"))(n=length(breaks)))

# get aggregate assessment
qc.stf.gseanno <- ifelse(st.agg$qc.stf.filt == "above", TRUE, FALSE)
ba.stf.gseanno <- ifelse(st.agg$ba.stf.filt == "above", TRUE, FALSE)
agg.stf <- ifelse(ba.stf.gseanno | qc.stf.gseanno, "above", "below")

# filter studies
gsefilt <- agg.stf == "above"
hm.baf.subset <- hm.baf[gsefilt,]
hm.qcf.subset <- hm.qcf[gsefilt,]
gsesize.subset <- gsesize[gsefilt]
ro <- rev(order(gsesize.subset))

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

hmlist <- h1 + h2 + hm.ranno

#pdf("hmcomp_subset_gse-stfreq.pdf", 7, 7)
#draw(hmlist, row_title = paste0("Study ID (N = ", length(gsesize.subset),")"), column_title = "Sub-threshold Frequency", 
#     heatmap_legend_side = "right", annotation_legend_side = "right")
#dev.off()