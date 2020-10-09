#!/usr/bin/env R

# Make figure 2A, BeadArray controls stacked barplot
# of sample performance outcomes.

library(ggplot2)
library(ggforce)
library(data.table)

library(recountmethylationManuscriptSupplement)
pkgname <- "recountmethylationManuscriptSupplement"
tables.dir <- system.file("extdata", "tables", package = pkgname)

#----------
# load data
#----------
fn <- "table-s2_qcmd-allgsm.csv"
qcmd <- fread(file.path(tables.dir, fn), sep = ",", data.table = FALSE)

which.bacol <- which(grepl("^ba.*", colnames(qcmd)))
bat <- qcmd[,c(1, which.bacol)]
colnames(bat) <- gsub("^ba\\.", "", colnames(bat))

#---------------
# make plot data
#---------------
tr <- bathresh(bat)
dfplot <- matrix(nrow = 0, ncol = 3)
for(c in colnames(tr)[2:18]){
  num.pass <- length(which(tr[,c] == "PASS"))
  num.fail <- length(which(tr[,c] == "FAIL"))
  dfplot = rbind(dfplot, matrix(c(num.pass, "ABOVE", c), nrow = 1))
  dfplot = rbind(dfplot, matrix(c(num.fail, "BELOW", c), nrow = 1))
}
dfplot = as.data.frame(dfplot, stringsAsFactors = F)
dfplot[,1] = as.numeric(dfplot[,1])
colnames(dfplot) = c("Num. Samples", "Thresh. Count", "Metric")
dfi = dfplot[dfplot[,2] == "BELOW",]
dfplot$Metric = factor(dfplot$Metric, levels = dfi[,3][rev(order(dfi[,1]))])

#----------
# make plot
#----------
fig2a <- ggplot(dfplot, aes(x = `Metric`, y = `Num. Samples`, fill = `Thresh. Count`)) +
  geom_bar(stat = 'identity') + scale_fill_manual(values = c("limegreen", "red")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_zoom(ylim = c(0, 2500), zoom.data = ifelse(a <= 10, NA, FALSE)) +
  xlab("BeadArray Metric") + ylab("Number of\nSamples")
