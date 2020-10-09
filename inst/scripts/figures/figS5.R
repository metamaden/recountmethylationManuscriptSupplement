#!/usr/bin/env R

# Make panel Fig 2C, study-level quality assessment outcomes

library(ComplexHeatmap)
library(circlize)

#----------
# load data
#----------
pkgname <- "recountmethylationManuscriptSupplement"
tables.dir <- system.file("extdata", "tables", package = pkgname)
tname <- "table-s5_qcfreq_gsewise.csv"
st.agg <- read.csv(file.path(tables.dir, tname), header = TRUE)

#-----------------
# get heatmap data
#-----------------
# heatmap data
hm.baf <- st.agg[,3:19]
hm.qcf <- st.agg[,20:25]
gsesize <- st.agg$ngsm
names(gsesize) <- st.agg$gseid
# reorder on size
ro <- rev(order(gsesize))
# get heatmap color palettes
breaks1 <- seq(min(hm.baf), max(hm.baf), 0.01)
hmcol1 <- colorRamp2(breaks1, colorRampPalette(c("blue","orange","red"))(n=length(breaks1)))
breaks2 <- seq(min(hm.qcf), max(hm.qcf), 0.01)
hmcol2 <- colorRamp2(breaks2, colorRampPalette(c("black","green","white"))(n=length(breaks2)))

# get hi stf ba metrics
thresh <- 0.6
ba.stf.gseanno <- c()
for(g in 1:nrow(hm.baf)){
  ba.stf.gseanno <- c(ba.stf.gseanno,
                     ifelse(length(which(hm.baf[g,] > thresh)) > 0, 
                            "above", "below"))
}

thresh <- 0.6
qc.stf.gseanno <- c()
which.col <- c(3, 4)
for(g in 1:nrow(hm.baf)){
  qc.stf.gseanno <- c(qc.stf.gseanno,
                     ifelse(length(which(hm.qcf[g, which.col] > thresh)) == 2, 
                            "above", "below"))
}

#------------------
# make plot objects
#------------------
h1.anno <- rowAnnotation(df = data.frame("ba.thresh" = ba.stf.gseanno[ro]), 
                        col = list("ba.thresh" = c("above" =  "red", "below" = "forestgreen")), 
                        width = unit(1, "cm"))

h2.anno <- rowAnnotation(df = data.frame("qc.thresh" = qc.stf.gseanno[ro]), 
                        col = list("qc.thresh" = c("above" =  "yellow", "below" = "purple")), 
                        width = unit(1, "cm"))

agg.stf <- ifelse(ba.stf.gseanno == "above" | qc.stf.gseanno == "above",
                 "above", "below")

h3.anno <- rowAnnotation(df = data.frame("agg.thresh" = agg.stf[ro]), 
                        col = list("agg.thresh" = c("above" =  "white", "below" = "black")), 
                        width = unit(1, "cm"))

h1 <- Heatmap(hm.baf[ro,], name = "BeadArray\nS.T. Freq.",
             col = hmcol1, cluster_rows = F, cluster_columns = T,
             show_row_dend = F, show_column_dend = F,
             show_row_names = F)

h2 <- Heatmap(hm.qcf[ro,], name = "Signal\nS.T. Freq.",
             col = hmcol2,
             cluster_rows = F, cluster_columns = T,
             show_row_dend = F, show_column_dend = F,
             show_row_names = F)


hm.ranno <- rowAnnotation(`Log2 Study Size\n(Num. GSM)` = anno_barplot(log2(gsesize[ro]), 
                                                                      width = unit(4, "cm")))

# objects as single plot list
hmlist <- h1.anno + h1 + h2.anno + h2 + h3.anno + hm.ranno

#--------------------
# make composite plot
#--------------------
row.title <- paste0("Study ID (N = ", length(gsesize),")")
col.title <- "Sub-threshold Frequency"
draw(hmlist, row_title = row.title, column_title = col.title, 
     heatmap_legend_side = "right", annotation_legend_side = "right")
