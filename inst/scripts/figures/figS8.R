#!/usr/bin/env

# Author: Sean Maden
# Summarize genome region mappings at probes with
# low variance and mean DNAm differences/ranges across
# 7 noncancer tissues (figS8).

library(ggplot2)
library(gridExtra)
library(minfiData)

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
cgaf <- cga[cga$Name %in% cgid$cgid,]
mainplotwidth = 8; mainplotheight = 8
theme.size = 20

# plot 1 -- island/gene region
bpdf <- as.data.frame(table(cgaf$annocat), stringsAsFactors = FALSE)
colnames(bpdf) = c("Genome\nRegion", "Count")
bpdf[,1] <- ifelse(bpdf[,1] == "interisland_intergenic", "opensea\n",
                   ifelse(bpdf[,1] == "interisland_intragenic", "opensea;\nintragenic",
                          ifelse(bpdf[,1] == "intraisland_intergenic", "island\n",
                                 "island;\nintragenic")))

figS8.plot1 <- ggplot(bpdf, aes(x = "", y = Count, fill = `Genome\nRegion`)) +
  geom_bar(stat = "identity") + theme_bw(base_size = theme.size) +
  theme(axis.text.x = element_text(angle = 90), 
        legend.position = "top") + 
  xlab("") + ylab("Number of Probes") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# plot 2 -- gene region
bpdf <- matrix(nrow = 0, ncol = 2)
bpdf <- as.data.frame(table(cgaf$gene.relation),
                     stringsAsFactors = FALSE)
bpdf[,1] <- ifelse(bpdf[,1] == "body", "body\n",
                   ifelse(bpdf[,1] == "intergenic", "intergenic\n",
                          ifelse(bpdf[,1] == "promoter", "promoter\n",
                                 "promoter;\nbody")))
colnames(bpdf) <- c("Gene\nRegion", "Count")
figS8.plot2 <- ggplot(bpdf, aes(x = "", y = Count, fill = `Gene\nRegion`)) +
  geom_bar(stat = "identity") + theme_bw(base_size = theme.size) +
  theme(axis.text.x = element_text(angle = 90), 
        legend.position = "top") + 
  xlab("") + ylab("") +
  guides(fill = guide_legend(nrow=2, byrow=TRUE))

# plot 3 -- island region
bpdf <- matrix(nrow = 0, ncol = 2)
cgaf$isl.bin <- ifelse(grepl("Shore", cgaf$Relation_to_Island), "shore\n",
                       ifelse(grepl("Shelf", cgaf$Relation_to_Island), "shelf\n",
                              paste0(tolower(cgaf$Relation_to_Island), "\n")))
bpdf <- as.data.frame(table(cgaf$isl.bin), stringsAsFactors = FALSE)
colnames(bpdf) <- c("Island\nRegion", "Count")
figS8.plot3 <- ggplot(bpdf, aes(x = "", y = Count, fill = `Island\nRegion`)) +
  geom_bar(stat = "identity") + theme_bw(base_size = theme.size) +
  theme(axis.text.x = element_text(angle = 90), 
        legend.position = "top") + 
  xlab("") + ylab("") +
  guides(fill = guide_legend(nrow=2,byrow=TRUE))

#-------------------------------
# make horizontal composite plot
#-------------------------------
#bpstack.name <- "bp-region_lowvar_nct7.pdf"
#pdf(bpstack.name, 15, 5.5)
#grid.arrange(figS8.plot1, figS8.plot2, figS8.plot3, ncol = 3)
#dev.off()