#!/usr/bin/env R

# Make Figure 1A panel plot of cumulative GEO samples (GSMs)
# by year and available sample type.

library(ggplot2)
library(data.table)

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
fnv <- c("gsmyeardata", "gsmidatyrdat")
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

#------------
# make ggplot
#------------
# make plot
dff$`Platform (Type)` <- dff$plab
lv <- c("HM450K (all)", "HM450K (idat)", "EPIC (all)", "EPIC (idat)", "HM27K (all)", "HM27K (idat)")
vv <- c("green", "forestgreen", "blue", "cyan", "gold", "brown")
sv <- c(17, 17, 16, 16, 18, 18)
lname <- "Platform (Type)"
fig1a <- ggplot(dff, aes(year, samples)) + theme_bw() +
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