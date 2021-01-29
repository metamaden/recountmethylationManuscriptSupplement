#!/usr/bin/env R

# Author: Sean Maden
# Get genome regions in high variance probes across 7 noncancer tissues, 
# then map with stacked barplots (fig4c).

library(minfiData)

#----------
# load data
#----------
pkgname <- "recountmethylationManuscriptSupplement"
tables.dir <- system.file("extdata", "tables", package = pkgname)
nct7.dir <- system.file("extdata", "nct7", package = pkgname)
ldfn <- "ltxfilt-dnamsummary_highvar_nct7.rda"
txfiltl <- get(load(file.path(nct7.dir, ldfn)))
txnames <- names(txfiltl)

# get minfi-formatted probe annotations
cga <- getAnnotation(get(data("MsetEx")))

#---------------------
# add modified cg anno
#---------------------
prom.stat <- grepl("TSS|5'", cga$UCSC_RefGene_Group)
body.stat <- grepl("Body|Exon|3'", cga$UCSC_RefGene_Group)
cga$gene.type <- ifelse(cga$UCSC_RefGene_Name=="", "intergenic",
                       ifelse(prom.stat & body.stat, 
                              "intragenic_promoter-body",
                              ifelse(prom.stat & !body.stat, 
                                     "intragenic_promoter", 
                                     "intragenic_body")))
cga$isl.type <- ifelse(cga$Relation_to_Island=="OpenSea", 
                       "interisland_opensea", 
                      ifelse(cga$Relation_to_Island=="Island", 
                             "intraisland_main", 
                             "intraisland_other"))
cga$type.composite <- paste0(cga$isl.type,";",cga$gene.type)
cga$annocat <- ifelse(cga$Relation_to_Island=="OpenSea" & 
                        cga$UCSC_RefGene_Name=="", "interisland_intergenic",
                     ifelse(!(cga$Relation_to_Island=="OpenSea") & 
                              cga$UCSC_RefGene_Name=="", 
                            "intraisland_intergenic",
                            ifelse(cga$Relation_to_Island=="OpenSea" & 
                                     !(cga$UCSC_RefGene_Name==""), 
                                   "interisland_intragenic",
                                   ifelse(!(cga$Relation_to_Island=="OpenSea" | 
                                              cga$UCSC_RefGene_Name==""), 
                                          "intraisland_intragenic", "NA"))))
prom.stat <- grepl("TSS|5'", cga$UCSC_RefGene_Group)
body.stat <- grepl("Body|Exon|3'", cga$UCSC_RefGene_Group)
cga$gene.relation = ifelse(cga$UCSC_RefGene_Name=="", "intergenic",
                           ifelse(prom.stat & body.stat, "promoter;body",
                                  ifelse(prom.stat & !body.stat, "promoter", "body")))

#------------------------
# barplots, mapping freqs
#------------------------
txnames <- names(txfiltl);mainplotwidth <- 8; mainplotheight <- 8
legendpos <- "top"

# top 1k mvp, tissue-specific
bpdf <- matrix(nrow = 0, ncol = 3);which.lltx <- 6
for(ti in 1:length(txnames)){
  txname <- txnames[ti]; mvpidt = txfiltl[[txname]]$cgid
  tdati <- as.data.frame(table(cga[rownames(cga) %in% mvpidt,]$annocat))
  tdati$txname <- txname; bpdf = rbind(bpdf, tdati)
}; bpdf <- as.data.frame(bpdf, stringsAsFactors = FALSE)
colnames(bpdf) <- c("Genomic Region", "Count", "Tissue")
bpdf$label <- as.character(bpdf$`Genomic Region`)
bpdf$Tissue <- factor(bpdf$Tissue, levels=txnames)

# get compact label from region var
bpdf[bpdf$label == "interisland_intergenic",]$label <- "Open sea,\nintergenic"
bpdf[bpdf$label == "interisland_intragenic",]$label <- "Open sea,\nintragenic"
bpdf[bpdf$label == "intraisland_intergenic",]$label <- "CpG island,\nintergenic"
bpdf[bpdf$label == "intraisland_intragenic",]$label <- "CpG island,\nintragenic"

# font sizes 
fs.axis.text <- fs.legend.text <- 10
fs.axis.title <- fs.legend.title <- 12

# new plot, compact legend format
p1 <- ggplot(bpdf, aes(x = Tissue, y = Count, fill = label)) +
  geom_bar(stat = "identity", colour = "black") + theme_bw() + 
  labs(fill = "Genome\nregion") + ylab("Number of probes") + xlab("") +
  theme(axis.text.x = element_text(angle = 90), legend.position = legendpos, 
        axis.text = element_text(size = fs.axis.text), 
        axis.title = element_text(size = fs.axis.title),
        legend.text = element_text(size = fs.legend.text),
        legend.title = element_text(size = fs.legend.title)) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
  

#--------------------------
# anno plots: gene relation
#--------------------------
bpdf <- matrix(nrow = 0, ncol = 3);which.lltx <- 6
for(ti in 1:length(txnames)){
  txname <- txnames[ti]; mvpidt = txfiltl[[txname]]$cgid
  tdati <- as.data.frame(table(cga[rownames(cga) %in% 
                                     mvpidt,]$gene.relation))
  tdati$txname <- txname; bpdf = rbind(bpdf, tdati)
};bpdf <- as.data.frame(bpdf, stringsAsFactors = FALSE)
colnames(bpdf) <- c("Gene Region", "Count", "Tissue")
bpdf$Tissue <- factor(bpdf$Tissue, levels=txnames)
bpdf$label <- as.character(bpdf$`Gene Region`)

# modify labels with exclusions
bpdf[bpdf$label == "promoter;body",]$label <- "Promoter and\nbody"
bpdf[bpdf$label == "body",]$label <- "Body only\n "
bpdf[bpdf$label == "promoter",]$label <- "Promoter\nonly"
bpdf[bpdf$label == "intergenic",]$label <- "Intergenic\n "

# font sizes 
fs.axis.text <- fs.legend.text <- 10
fs.axis.title <- fs.legend.title <- 12

# new plot with compact legend
p2 <- ggplot(bpdf, aes(x = Tissue, y = Count, fill = label)) +
  geom_bar(stat = "identity", colour = "black") + theme_bw() + 
  labs(fill = "Gene\nregion") + xlab("Tissue") +
  theme(axis.text.x = element_text(angle = 90), legend.position = legendpos,
        axis.title.y = element_blank(), 
        axis.text = element_text(size = fs.axis.text), 
        axis.title = element_text(size = fs.axis.title), 
        legend.text = element_text(size = fs.legend.text),
        legend.title = element_text(size = fs.legend.title)) + 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

#-------------------------
# anno plots: isl relation
#-------------------------
# top 1k mvp, tissue-specific
bpdf <- matrix(nrow = 0, ncol = 3);which.lltx <- 6
for(ti in 1:length(txnames)){
  txname <- txnames[ti]; mvpidt = txfiltl[[txname]]$cgid
  tdati <- as.data.frame(table(cga[rownames(cga) %in% 
                                     mvpidt,]$Relation_to_Island))
  tdf <- tdati[grepl("Island|OpenSea", tdati[,1]),]
  shore.sum <- sum(tdati[grepl("Shore", tdati[,1]),2])
  shelf.sum <- sum(tdati[grepl("Shelf", tdati[,1]),2])
  tdi1 <- matrix(c("shore", shore.sum), nrow = 1)
  tdi2 <- matrix(c("shelf", shelf.sum), nrow = 1)
  colnames(tdf) <- colnames(tdi1) <- colnames(tdi2) <- colnames(tdati)
  tdati <- rbind(tdf, rbind(tdi1, tdi2))
  tdati$txname <- txname; bpdf = rbind(bpdf, tdati)
};bpdf <- as.data.frame(bpdf, stringsAsFactors = FALSE)
colnames(bpdf) <- c("Island Region", "Count", "Tissue")
bpdf$label <- as.character(bpdf$`Island Region`)
bpdf$Tissue <- factor(bpdf$Tissue, levels=txnames)
bpdf$Count <- as.numeric(bpdf$Count)

# new compact label from region var
bpdf[bpdf$label == "Island",]$label <- "CpG island\n"
bpdf[bpdf$label == "OpenSea",]$label <- "Open sea\n"
bpdf[bpdf$label == "shelf",]$label <- "CpG island\nshelf"
bpdf[bpdf$label == "shore",]$label <- "CpG island\nshore"

# font sizes 
fs.axis.text <- fs.legend.text <- 10
fs.axis.title <- fs.legend.title <- 12

# new plot with compact legend
p3 <- ggplot(bpdf, aes(x = Tissue, y = Count, fill = label)) +
  geom_bar(stat = "identity", colour = "black") + theme_bw() + 
  labs(fill = "Island\nregion") + xlab("") +
  theme(axis.text.x = element_text(angle = 90), legend.position = legendpos,
        axis.title.y = element_blank(), axis.text = element_text(size = fs.axis.text), 
        axis.title = element_text(size = fs.axis.title), 
        legend.text = element_text(size = fs.legend.text),
        legend.title = element_text(size = fs.legend.title)) + 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

#---------------------------------------
# composite plot, horizontal arrangement
#---------------------------------------
# print for manuscript
#pdf("fig4c_bp-region_hivar_nct7.pdf", 10.5, 3.1)
#print(grid.arrange(p1, p2, p3, ncol = 3));dev.off()

# print for vignette
print(grid.arrange(p1, p2, p3, ncol = 3))