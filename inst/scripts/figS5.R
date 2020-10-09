#!/usr/bin/env R

# Make panel Fig 2C, study-level quality assessment outcomes

library(ComplexHeatmap)
library(circlize)

#----------
# load data
#----------

library(recountmethylationManuscriptSupplement)
pkgname <- "recountmethylationManuscriptSupplement"
tables.dir <- system.file("extdata", "tables", package = pkgname)

#----------
# make plot
#----------
# reorder on size
identical(rownames(hm.baf), rownames(hm.qcf))
ro = rev(order(gsesize))

hmdat = hm.baf
breaks=seq(min(hmdat), max(hmdat), 0.01)
hmcol1 = colorRamp2(breaks, colorRampPalette(c("blue","orange","red"))(n=length(breaks)))

hmdat = hm.qcf
breaks=seq(min(hmdat), max(hmdat), 0.01)
hmcol2 = colorRamp2(breaks, colorRampPalette(c("black","green","white"))(n=length(breaks)))

# save the aggregate dataset


# get hi stf ba metrics
thresh = 0.6
ba.stf.gseanno = c()
for(g in 1:nrow(hm.baf)){
  ba.stf.gseanno = c(ba.stf.gseanno,
                     ifelse(length(which(hm.baf[g,] > thresh)) > 0, 
                            "above", "below"))
}

thresh = 0.6
qc.stf.gseanno = c()
which.col = c(3, 4)
for(g in 1:nrow(hm.baf)){
  qc.stf.gseanno = c(qc.stf.gseanno,
                     ifelse(length(which(hm.qcf[g, which.col] > thresh)) == 2, 
                            "above", "below"))
}

h1.anno = rowAnnotation(df = data.frame("ba.thresh" = ba.stf.gseanno[ro]), 
                        col = list("ba.thresh" = c("above" =  "red", "below" = "forestgreen")), 
                        width = unit(1, "cm"))

h2.anno = rowAnnotation(df = data.frame("qc.thresh" = qc.stf.gseanno[ro]), 
                        col = list("qc.thresh" = c("above" =  "yellow", "below" = "purple")), 
                        width = unit(1, "cm"))

agg.stf = ifelse(ba.stf.gseanno == "above" | qc.stf.gseanno == "above",
                 "above", "below")

h3.anno = rowAnnotation(df = data.frame("agg.thresh" = agg.stf[ro]), 
                        col = list("agg.thresh" = c("above" =  "white", "below" = "black")), 
                        width = unit(1, "cm"))

h1 = Heatmap(hm.baf[ro,], name = "BeadArray\nS.T. Freq.",
             col = hmcol1, cluster_rows = F, cluster_columns = T,
             show_row_dend = F, show_column_dend = F,
             show_row_names = F)

h2 = Heatmap(hm.qcf[ro,], name = "Signal\nS.T. Freq.",
             col = hmcol2,
             cluster_rows = F, cluster_columns = T,
             show_row_dend = F, show_column_dend = F,
             show_row_names = F)


hm.ranno = rowAnnotation(`Log2 Study Size\n(Num. GSM)` = anno_barplot(log2(gsesize[ro]), 
                                                                      width = unit(4, "cm")))

hmlist = h1.anno + h1 + h2.anno + h2 + h3.anno + hm.ranno


pdf("hmcomp_gse-stfreq.pdf", 8, 8)
draw(hmlist, row_title = paste0("Study ID (N = ", length(gsesize),")"), column_title = "Sub-threshold Frequency", 
     heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()