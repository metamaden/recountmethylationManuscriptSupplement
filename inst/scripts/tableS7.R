#!/usr/bin/env

# Summarize probe and genome region mappings at probes
# with low variances and mean ranges across 7 noncancer 
# tissues (Table S7).

#----------
# load data
#----------
# get minfi-formatted probe annotations
cga <- getAnnotation(get(data("MsetEx")))

# get probes
ptidf <- get(load(ptidfn))
st <- data.frame(cgid = as.character(ptidf), 
                 stringsAsFactors = FALSE)

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

#---------------
# make new table
#---------------
# get anno and match data
cgf = cga[cga$Name %in% ptidf,]
cgf = cgf[order(match(cgf$Name, st$cgid)),]
identical(cgf$Name, st$cgid)
sta = cbind(st, cgf[,c(1, 2, 3, 9, 14, 15, 20, 21, 22)])

# save
#sta.name <- "supptable_recurr-lowvar-cg-anno_7nct"
#save(sta, file = paste0(sta.name, ".rda"))
#write.csv(sta, file = paste0(sta.name, ".csv"), row.names = FALSE)