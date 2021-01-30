#!/usr/bin/env R

# Author: Sean Maden
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
# format colnames 
cv <- colnames(st.agg); which.cv <- grepl("meth|unmeth", cv)
colnames(st.agg)[which.cv] <- gsub("_", " ", 
                                   gsub("\\.", "<", 
                                        colnames(st.agg)[which.cv]))
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
hmcol2 <- colorRamp2(breaks2, colorRampPalette(c("black","forestgreen","green"))(n=length(breaks2)))

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
h1.legend.title <- list(title=expression(italic("BeadArray f"["st"])))
h2.legend.title <- list(title=expression(italic("Signal f"["st"])))

h1.anno <- rowAnnotation(df = data.frame("BeadArray_thresh." = ba.stf.gseanno[ro]), 
                        col = list("BeadArray_thresh." = c("above" =  "red", 
                                                           "below" = "forestgreen")), 
                        width = unit(1, "cm"))

h2.anno <- rowAnnotation(df = data.frame("Signal_thresh." = qc.stf.gseanno[ro]), 
                        col = list("Signal_thresh." = c("above" =  "yellow", 
                                                        "below" = "purple")), 
                        width = unit(1, "cm"))

agg.stf <- ifelse(ba.stf.gseanno == "above" | qc.stf.gseanno == "above",
                 "above", "below")

h3.anno <- rowAnnotation(df = data.frame("Aggregate_thresh." = agg.stf[ro]), 
                        col = list("Aggregate_thresh." = c("above" =  "white", 
                                                           "below" = "black")), 
                        width = unit(1, "cm"))

h1 <- Heatmap(hm.baf[ro,], heatmap_legend_param = h1.legend.title,
             col = hmcol1, cluster_rows = F, cluster_columns = T,
             show_row_dend = F, show_column_dend = F,
             show_row_names = F)

h2 <- Heatmap(hm.qcf[ro,], heatmap_legend_param = h2.legend.title,
             col = hmcol2, cluster_rows = F, cluster_columns = T,
             show_row_dend = F, show_column_dend = F,
             show_row_names = F)


hm.ranno <- rowAnnotation(`Log2 study size\n(number of samples)` = anno_barplot(log2(gsesize[ro]), 
                                                                      width = unit(4, "cm")))

# objects as single plot list
hmlist <- h1.anno + h1 + h2.anno + h2 +  hm.ranno + h3.anno

#--------------------
# make composite plot
#--------------------
row.title <- paste0(paste(rep(" ", 30), collapse = ""), 
                    "Study (",length(gsesize)," studies)"); col.title <- ""

# print manuscript plot
#pdf("sfig5_hmcomp_gse-stfreq.pdf", 8.8, 4)
#draw(hmlist, row_title = row.title, column_title = col.title, 
#     heatmap_legend_side = "right", annotation_legend_side = "right")
#dev.off()

# draw vignette plot
draw(hmlist, row_title = row.title, column_title = col.title, 
     heatmap_legend_side = "right", annotation_legend_side = "right")
