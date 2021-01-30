#!/usr/bin/env R

# Author: Sean Maden
# Make Figure 1A panel plot of cumulative GEO samples (GSMs) by year 
# and available sample type.

library(ggplot2)
library(data.table)
library(ggpubr)
library(cowplot)

pkgname <- "recountmethylationManuscriptSupplement"
tables.dir <- system.file("extdata", "geodata", package = pkgname)

#--------------
# platform info
#--------------
lid <- list(c('GPL13534','HM450K','2011','May 2011'),
            c('GPL21145', 'EPIC','2015','Nov. 2015'),
            c('GPL8490', 'HM27K','2009','Apr. 2009')
)
keydf <- as.data.frame(do.call(rbind, lid))
colnames(keydf) <- c("platform_id", "platform_name", 
                     "release_year", "release_month")

#----------------------------
# load and format data tables
#----------------------------
# specify tables
fnv <- c("gseyeardata", "gseidatyrdat")
# load sample data
dm <- do.call(rbind, lapply(fnv, function(x){
  return(fread(file.path(tables.dir, x), sep = " ", data.table = FALSE))
}))
df <- as.data.frame(dm, stringsAsFactors = FALSE)
colnames(df) <- c("platformid","year","samples")
# filter platforms of interest
which.platforms <- c("GPL8490", "GPL13534", "GPL21145")
df <- df[df$platformid %in% which.platforms,]
# append sample types
stype <- c(rep(" (all)", 60), rep(" (idat)", 60))
df$plab <- paste0(df$platformid, stype)
# get formatted columns/labels
df$year <- as.numeric(df$year); df$pname <- "NA"
dff <- matrix(nrow = 0, ncol = 5)
for(pid in keydf[,1]){
  ki <- keydf[keydf[,1] == pid,]
  release.yr <- as.numeric(ki[3])
  dpid <- df[df$platformid == pid,]
  dpid <- dpid[dpid$year >= release.yr,]
  dpid$pname <- paste0(gsub(" .*", "", ki[2]),
                       " ", gsub(".* ", "", dpid$plab))
  dpid <- dpid[order(dpid$year),]
  for(u in unique(dpid$plab)){
    dpp <- dpid[dpid$plab == u,]
    for(r in seq(nrow(dpp))){
      dfi <- dpp[r,] 
      dfi[3] <- sum(dpp[seq(r),3])
      dff <- rbind(dff, dfi)
    }
  }
}
dff <- as.data.frame(dff, stringsAsFactors = FALSE)
dff$platform <- gsub(" .*", "", dff$pname)

#---------------------
# make manuscript plot
#---------------------
# make manuscript plot
dff$`Platform (Type)` <- dff$plab
lv <- c("HM450K (all)", "HM450K (idat)", "EPIC (all)", "EPIC (idat)", "HM27K (all)", "HM27K (idat)")
vv <- c("green", "forestgreen", "blue", "cyan", "gold", "brown")
sv <- c(17, 17, 16, 16, 18, 18)
lname <- ""; lsize <- 1; ptsize <- 4.5
yearv <- c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018)

vv2 <- c("HM27K (all)" = "brown", "HM27K (idat)" = "gold",
         "HM450K (all)" = "green", "HM450K (idat)" = "forestgreen",
         "EPIC (all)" = "cyan", "EPIC (idat)" = "blue")
sv2 <- c("HM27K (all)" = 18, "HM27K (idat)" = 18,
         "HM450K (all)" = 16, "HM450K (idat)" = 16,
         "EPIC (all)" = 17, "EPIC (idat)" = 17)

dff$pname <- factor(dff$pname, levels = c("HM27K (all)","HM27K (idat)",
                                          "HM450K (all)", "HM450K (idat)",
                                          "EPIC (all)", "EPIC (idat)"))
dff$color <- dff$shape <- "NA"
for(x in names(vv2)){
  dff[dff$pname == x,]$color <- vv2[x];dff[dff$pname == x,]$shape <- sv2[x]}

figS1 <- ggplot(dff, aes(x = year, y = samples)) + theme_bw() +
  xlab("Year") + ylab("Cumulative studies") + xlim(2009, 2018) +
  geom_line(data = dff[dff$pname == "HM27K (all)",], 
            aes(year, samples), color = "brown", size = lsize) +
  geom_point(data = dff[dff$pname == "HM27K (all)",], aes(year, samples), 
             color = "brown", size = ptsize, shape = 18) +
  geom_line(data = dff[dff$pname == "HM27K (idat)",], 
            aes(year, samples), color = "gold", size = lsize) +
  geom_point(data = dff[dff$pname == "HM27K (idat)",], aes(year, samples), 
             color = "gold", size = ptsize, shape = 18) +
  geom_line(data = dff[dff$pname == "HM450K (all)",], 
            aes(year, samples), color = "green", size = lsize) +
  geom_point(data = dff[dff$pname == "HM450K (all)",], aes(year, samples), 
             color = "green", size = ptsize, shape = 16) +
  geom_line(data = dff[dff$pname == "HM450K (idat)",], 
            aes(year, samples), color = "forestgreen", size = lsize) +
  geom_point(data = dff[dff$pname == "HM450K (idat)",], aes(year, samples), 
             color = "forestgreen", size = ptsize, shape = 16) +
  geom_line(data = dff[dff$pname == "EPIC (all)",], 
            aes(year, samples), color = "cyan", size = lsize) +
  geom_point(data = dff[dff$pname == "EPIC (all)",], aes(year, samples), 
             color = "cyan", size = ptsize, shape = 17) +
  geom_line(data = dff[dff$pname == "EPIC (idat)",], 
            aes(year, samples), color = "blue", size = lsize) +
  geom_point(data = dff[dff$pname == "EPIC (idat)",], aes(year, samples), 
             color = "blue", size = ptsize, shape = 17) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
  scale_x_continuous(labels = yearv,
                     breaks = yearv)

# get the plot legend
figS1.lform <- ggplot(dff, aes(x = year, y = samples)) + theme_bw() +
  xlab("Year") + ylab("Cumulative samples") + xlim(2009, 2018) +
  geom_point(aes(colour = pname, shape = pname, size = 4)) +
  geom_line(aes(colour = pname, size = 2)) +
  scale_colour_manual(values = vv2) + scale_shape_manual(values = sv2) +
  guides(colour = FALSE,
         shape = guide_legend(override.aes = list(shape = sv2, 
                                                  colour = vv2, size = 4)), 
         size = FALSE) +
  theme(legend.title = element_blank(), legend.position = c(0.2, 0.6))

figS1.legend <- ggpubr::get_legend(figS1.lform)

# print manuscript plot
#pdf("sfig1_geo-gse-yr.pdf", 5, 3.1)
#print(cowplot::ggdraw(figS1) + 
#        cowplot::draw_plot(figS1.legend, .17, .38, .5, .5))
#dev.off()

#-------------------
# make vignette plot
#-------------------
# make plot
dff$`Platform (Type)` <- dff$plab
lv <- c("HM450K (all)", "HM450K (idat)", "EPIC (all)", "EPIC (idat)", "HM27K (all)", "HM27K (idat)")
vv <- c("green", "forestgreen", "blue", "cyan", "gold", "brown")
sv <- c(17, 17, 16, 16, 18, 18)
lname <- "Platform (Type)"
figS1 <- ggplot(dff, aes(year, samples)) + theme_bw() +
  xlab("Year") + ylab("Cumulative Samples") + 
  geom_point(aes(colour = `Platform (Type)`, shape = platform)) + 
  guides(shape = FALSE) +
  scale_shape_manual(values = sv) + 
  scale_colour_manual(values = vv) +
  guides(colour = guide_legend(override.aes = list(shape = sv, colour = vv))) +
  geom_line(data = dff[dff$pname == "HM27K (idat)",], 
            aes(year, samples), color = "brown") +
  geom_line(data = dff[dff$pname == "HM27K (all)",], 
            aes(year, samples), color = "gold") +
  geom_line(data = dff[dff$pname == "HM450K (all)",], 
            aes(year, samples), color = "green") +
  geom_line(data = dff[dff$pname == "HM450K (idat)",], 
            aes(year, samples), color = "forestgreen") +
  geom_line(data = dff[dff$pname == "EPIC (all)",], 
            aes(year, samples), color = "blue") +
  geom_line(data = dff[dff$pname == "EPIC (idat)",], 
            aes(year, samples), color = "cyan") +
  theme(legend.position = c(0.2, 0.6))