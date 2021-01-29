#!/usr/bin/env R

# Script to recreate Fig 2B storage condition quality signals plot.

library(ggplot2); library(gridExtra); library(recountmethylation)

#----------
# load data
#----------
sf <- system.file(file.path("extdata", "data_analyses"), 
                  package = "recountmethylation")
load(file.path(sf, "data_analyses.RData"))

#----------------------
# format plot variables
#----------------------
ds[ds$storage == "ffpe",]$storage <- "FFPE"
ds[ds$storage == "frozen",]$storage <- "Fresh frozen"
colnames(ds)[3] <- "Storage\ncondition"

#----------------------
# make individual plots
#----------------------
fig2b.gp <- ggplot(ds, aes(x = meth.l2med, y = unmeth.l2med, color = `Storage\ncondition`)) + 
  geom_point(alpha = 0.35, cex = 3) + theme_bw() +
  scale_color_manual(values = c("FFPE" = "orange", "Fresh frozen" = "purple")) +
  xlab("Methylated signal (log2 scale)") + ylab("Unmethylated signal (log2 scale)")

fig2b.ellipse <- ggplot(ds, aes(x = meth.l2med, y = unmeth.l2med, color = `Storage\ncondition`)) + 
  stat_ellipse(level = 0.95) + theme_bw() +
  scale_color_manual(values = c("FFPE" = "orange", "Fresh frozen" = "purple")) +
  xlab("Methylated signal (log2 scale)") + ylab("Unmethylated signal (log2 scale)")

#------------------------------
# make formatted composite plot
#------------------------------
# font sizes
fs.axis.text <- fs.legend.text <- 11
fs.axis.title <- fs.legend.title <- 13

fig2b.gp <- ggplot(ds, aes(x = meth.l2med, y = unmeth.l2med, color = `Storage\ncondition`)) + 
  geom_point(alpha = 0.35, cex = 2) + theme_bw() + xlab("") + ylab("") +
  scale_color_manual(values = c("FFPE" = "orange", "Fresh frozen" = "purple")) +
  theme(legend.position = "none", plot.margin=unit(c(0.1, 0, -0.4, 0), "cm"),
        axis.text = element_text(size = fs.axis.text), 
        axis.title = element_text(size = fs.axis.title),
        legend.text = element_text(size = fs.legend.text),
        legend.title = element_text(size = fs.legend.title)) + 
  xlim(7, 14) + ylim(7, 14) + ylab("Unmethylated signal\n(log2 scale)")

fig2b.ellipse <- ggplot(ds, aes(x = meth.l2med, y = unmeth.l2med, color = `Storage\ncondition`)) + 
  stat_ellipse(level = 0.95) + theme_bw() + xlab("") + ylab("") +
  scale_color_manual(values = c("FFPE" = "orange", "Fresh frozen" = "purple")) +
  xlim(7, 14) + ylim(7, 14) + 
  theme(plot.margin=unit(c(0.1, 0, -0.4, 0), "cm"),
        axis.text = element_text(size = fs.axis.text), 
        axis.title = element_text(size = fs.axis.title),
        legend.text = element_text(size = fs.legend.text),
        legend.title = element_text(size = fs.legend.title))

lm <- matrix(c(rep(1,5),rep(2,7)),nrow = 1)

# print for manuscript
#pdf("fig2b_ggcomp_storage.pdf", 7.2, 3)
#print(grid.arrange(fig2b.gp, fig2b.ellipse, layout_matrix = lm,
#                   bottom = "Methylated signal\n(log2 scale)")); dev.off()

# print for vignette
grid.arrange(fig2b.gp, fig2b.ellipse, layout_matrix = lm,
             bottom = "Methylated signal\n(log2 scale)")