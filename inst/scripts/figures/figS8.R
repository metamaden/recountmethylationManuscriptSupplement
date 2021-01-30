#!/usr/bin/env

# Author: Sean Maden
# Summarize genome region mappings at probes with
# low variance and mean DNAm differences/ranges across
# 7 noncancer tissues (figS8).

library(ggplot2); library(gridExtra); library(minfiData)

#----------
# load data
#----------
pkgname <- "recountmethylationManuscriptSupplement"
tables.dir <- system.file("extdata", "tables", package = pkgname)
nct7.dir <- system.file("extdata", "nct7", package = pkgname)

# load data objects
lmvp.name = "lmvp_bin-abs-unique-lowvar_7nct.rda"
lmvp <- get(load(file.path(nct7.dir, lmvp.name)))
txnames <- names(lmvp) # tissues

# load the supplemental table
tname <- "table-s7_pantissue-hypovar_7nct.csv"
cgid <- read.csv(file.path(tables.dir, tname), header = TRUE)

# get minfi-formatted probe annotations
cga <- getAnnotation(get(data("MsetEx")))

#---------------------
# add modified cg anno
#---------------------
prom.stat = grepl("TSS|5'", cga$UCSC_RefGene_Group)
body.stat = grepl("Body|Exon|3'", cga$UCSC_RefGene_Group)
cga$gene.type = ifelse(cga$UCSC_RefGene_Name=="", "intergenic",
                       ifelse(prom.stat & body.stat, "intragenic_promoter-body",
                              ifelse(prom.stat & !body.stat, "intragenic_promoter", 
                                     "intragenic_body")))
cga$isl.type = ifelse(cga$Relation_to_Island=="OpenSea", "interisland_opensea", 
                      ifelse(cga$Relation_to_Island=="Island", "intraisland_main", 
                             "intraisland_other"))
table(cga$isl.type, cga$gene.type)
cga$type.composite = paste0(cga$isl.type,";",cga$gene.type)

cga$annocat = ifelse(cga$Relation_to_Island=="OpenSea" & cga$UCSC_RefGene_Name=="", 
                     "interisland_intergenic",
                     ifelse(!(cga$Relation_to_Island=="OpenSea") & cga$UCSC_RefGene_Name=="", 
                            "intraisland_intergenic",
                            ifelse(cga$Relation_to_Island=="OpenSea" & !(cga$UCSC_RefGene_Name==""), 
                                   "interisland_intragenic",
                                   ifelse(!(cga$Relation_to_Island=="OpenSea" | cga$UCSC_RefGene_Name==""), 
                                          "intraisland_intragenic", "NA"))))
prom.stat = grepl("TSS|5'", cga$UCSC_RefGene_Group)
body.stat = grepl("Body|Exon|3'", cga$UCSC_RefGene_Group)
cga$gene.relation = ifelse(cga$UCSC_RefGene_Name=="", "intergenic",
                           ifelse(prom.stat & body.stat, "promoter;body",
                                  ifelse(prom.stat & !body.stat, "promoter", "body")))

#---------------------------
# get region barplot objects
#---------------------------
cgaf <- cga[cga$Name %in% cgid$cgid,]; mainplotwidth = 8; mainplotheight = 8
theme.size = 20

# plot 1 -- island/gene region
bpdf <- as.data.frame(table(cgaf$annocat), stringsAsFactors = FALSE)
colnames(bpdf) = c("Genome\nregion", "Count"); bpdf$label <- bpdf[,1]
# get compact label from region var
bpdf[bpdf$label == "interisland_intergenic",]$label <- "Open sea,\nintergenic"
bpdf[bpdf$label == "interisland_intragenic",]$label <- "Open sea,\nintragenic"
bpdf[bpdf$label == "intraisland_intergenic",]$label <- "CpG island,\nintergenic"
bpdf[bpdf$label == "intraisland_intragenic",]$label <- "CpG island,\nintragenic"

figS8.plot1 <- ggplot(bpdf, aes(x = "", y = Count, 
                                fill = label)) +
  geom_bar(stat = "identity", colour = "black") + 
  theme_bw(base_size = theme.size) +
  theme(axis.text.x = element_text(angle = 90), 
        legend.position = "top") + xlab("") + ylab("Number of probes") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(fill = "Genome\nregion")

# plot 2 -- gene region
bpdf <- matrix(nrow = 0, ncol = 2)
bpdf <- as.data.frame(table(cgaf$gene.relation),
                     stringsAsFactors = FALSE);bpdf$label <- bpdf[,1]

bpdf[bpdf$label == "promoter;body",]$label <- "Promoter and\nbody"
bpdf[bpdf$label == "body",]$label <- "Body only\n "
bpdf[bpdf$label == "promoter",]$label <- "Promoter\nonly"
bpdf[bpdf$label == "intergenic",]$label <- "Intergenic\n "
colnames(bpdf) <- c("Gene\nregion", "Count", "label")

figS8.plot2 <- ggplot(bpdf, aes(x = "", y = Count, 
                                fill = label)) +
  geom_bar(stat = "identity", colour = "black") + 
  theme_bw(base_size = theme.size) +
  theme(axis.text.x = element_text(angle = 90), 
        legend.position = "top") + 
  xlab("") + ylab("") + labs(fill = "Gene\nregion") +
  guides(fill = guide_legend(nrow=2, byrow=TRUE))
  

# plot 3 -- island region
bpdf <- matrix(nrow = 0, ncol = 2)
cgaf$isl.bin <- ifelse(grepl("Shore", cgaf$Relation_to_Island), "shore",
                       ifelse(grepl("Shelf", cgaf$Relation_to_Island), "shelf",
                              paste0(tolower(cgaf$Relation_to_Island), "")))
bpdf <- as.data.frame(table(cgaf$isl.bin), stringsAsFactors = FALSE)
bpdf$label <- bpdf[,1];bpdf[bpdf$label == "island",]$label <- "CpG island\n"
bpdf[bpdf$label == "opensea",]$label <- "Open sea\n"
bpdf[bpdf$label == "shelf",]$label <- "CpG island\nshelf"
bpdf[bpdf$label == "shore",]$label <- "CpG island\nshore"
colnames(bpdf) <- c("Island\nregion", "Count", "label")

figS8.plot3 <- ggplot(bpdf, aes(x = "", y = Count, 
                                fill = label)) +
  geom_bar(stat = "identity", colour = "black") + 
  theme_bw(base_size = theme.size) +
  theme(axis.text.x = element_text(angle = 90), 
        legend.position = "top") + 
  xlab("") + ylab("") + labs(fill = "Island\nregion") +
  guides(fill = guide_legend(nrow=2,byrow=TRUE))
  
#-------------------------------
# make horizontal composite plot
#-------------------------------
#bpstack.name <- "sfig8_bp-region_lowvar_nct7.pdf"
#pdf(bpstack.name, 16, 5.5)
print(grid.arrange(figS8.plot1, figS8.plot2, figS8.plot3, ncol = 3))
#dev.off()