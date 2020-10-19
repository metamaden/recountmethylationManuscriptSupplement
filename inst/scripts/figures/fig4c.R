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
                       ifelse(prom.stat & body.stat, "intragenic_promoter-body",
                              ifelse(prom.stat & !body.stat, "intragenic_promoter", 
                                     "intragenic_body")))
cga$isl.type <- ifelse(cga$Relation_to_Island=="OpenSea", "interisland_opensea", 
                      ifelse(cga$Relation_to_Island=="Island", "intraisland_main", 
                             "intraisland_other"))
cga$type.composite <- paste0(cga$isl.type,";",cga$gene.type)
cga$annocat <- ifelse(cga$Relation_to_Island=="OpenSea" & 
                        cga$UCSC_RefGene_Name=="", "interisland_intergenic",
                     ifelse(!(cga$Relation_to_Island=="OpenSea") & 
                              cga$UCSC_RefGene_Name=="", "intraisland_intergenic",
                            ifelse(cga$Relation_to_Island=="OpenSea" & 
                                     !(cga$UCSC_RefGene_Name==""), "interisland_intragenic",
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
txnames <- names(txfiltl)
mainplotwidth <- 8; mainplotheight <- 8
legendpos <- "top"

# top 1k mvp, tissue-specific
bpdf <- matrix(nrow = 0, ncol = 3)
which.lltx <- 6
for(ti in 1:length(txnames)){
  txname <- txnames[ti]; mvpidt = txfiltl[[txname]]$cgid
  tdati <- as.data.frame(table(cga[rownames(cga) %in% mvpidt,]$annocat))
  tdati$txname <- txname; bpdf = rbind(bpdf, tdati)
}
colnames(bpdf) <- c("Genomic Region", "Count", "Tissue")
bpdf$Tissue <- factor(bpdf$Tissue, levels=txnames)

# get compact label from region var
bpdf$label <- "NA"
bpdf[bpdf$`Genomic Region` == "interisland_intergenic",]$label <- "opensea"
bpdf[bpdf$`Genomic Region` == "interisland_intragenic",]$label <- "opensea;intragenic"
bpdf[bpdf$`Genomic Region` == "intraisland_intergenic",]$label <- "island"
bpdf[bpdf$`Genomic Region` == "intraisland_intragenic",]$label <- "island;intragenic"

# new plot, compact legend format
p1 <- ggplot(bpdf, aes(x = Tissue, y = Count, fill = label)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.title = element_blank(),
        legend.position = legendpos,
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

#--------------------------
# anno plots: gene relation
#--------------------------
bpdf <- matrix(nrow = 0, ncol = 3)
which.lltx <- 6
for(ti in 1:length(txnames)){
  txname <- txnames[ti]; mvpidt = txfiltl[[txname]]$cgid
  tdati <- as.data.frame(table(cga[rownames(cga) %in% 
                                     mvpidt,]$gene.relation))
  tdati$txname <- txname; bpdf = rbind(bpdf, tdati)
}
colnames(bpdf) <- c("Gene Region", "Count", "Tissue")
bpdf$Tissue <- factor(bpdf$Tissue, levels=txnames)

# new compact label from region var
bpdf$label <- bpdf$`Gene Region`

# new plot with compact legend
p2 <- ggplot(bpdf, aes(x = Tissue, y = Count, fill = label)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = legendpos,
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) + 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

#-------------------------
# anno plots: isl relation
#-------------------------
# top 1k mvp, tissue-specific
bpdf <- matrix(nrow = 0, ncol = 3)
which.lltx <- 6
for(ti in 1:length(txnames)){
  txname <- txnames[ti]; mvpidt = txfiltl[[txname]]$cgid
  tdati <- as.data.frame(table(cga[rownames(cga) %in% 
                                     mvpidt,]$Relation_to_Island))
  tdati$txname <- txname; bpdf = rbind(bpdf, tdati)
}
colnames(bpdf) <- c("Island Region", "Count", "Tissue")
bpdf$Tissue <- factor(bpdf$Tissue, levels=txnames)

# new compact label from region var
bpdf$label <- tolower(bpdf$`Island Region`)
bpdf[bpdf$`Island Region` %in% c("N_Shelf", "S_Shelf"),]$label <- "shelf"
bpdf[bpdf$`Island Region` %in% c("N_Shore", "S_Shore"),]$label <- "shore"

# new plot with compact legend
p3 <- ggplot(bpdf, aes(x = Tissue, y = Count, fill = label)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = legendpos,
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) + 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

#---------------------------------------
# composite plot, horizontal arrangement
#---------------------------------------
grid.arrange(p1, p2, p3, ncol = 3, bottom = "Tissue", 
             left = "Number of\nProbes")
