#!/usr/bin/env R

library(ggplot2)
library(sp)
library(data.table)
library(gridExtra)
library(minfiData)

# This script describes variance analysis on the ANOVA-filtered
# probe data from 7 tissues (). First, run the script `nct_anova_bc.R`
# from a directory containing the `cgall` and `anova_results` folders
# with tissue-level ANOVA results and cg Beta-values summaries, which 
# generates `data.table.filt` files for each tissue. Next, run this 
# script to identify probes with cross-tissue low variances and 
# within-tissue high variances.

cga <- getAnnotation(get(data("MsetEx")))

lfilt.name = "lfilt_cgtables_nct.rda" # anova/chr filtered cg lists by tissue
lmvp.bin.name = "lmvp_01bin-quantfilt_nct.rda" # binned quantile mvp lists
lmvp.abs.name = "lmvp-all_5int_abs.rda" # absolute/global quantile mvp lists
lmvp.name = "lmvp_bin-abs-unique_7tissue.rda" # hypervar cgids by tissue, 2 methods
txmvpl.name = "ldat-txmvp_abs-bin_top1-5-10_7nct_" # fn stem for mvp lists by tissues
txfiltl.name = "listdat-top1kmvp-tx_abs-bin_7nct.rda" # filtered hypervar cgids data, 2 methods
sta.name = "supptable_1khypervar_7nct" # hypervar cg data supp table
st.name = "supptable_pantissue-hypovar_7nct" # hypovar cg data supp table
vp.name = "vp-1kmvp-tx_bin-abs_7nct.pdf" # violin plots fn
bpstack.name = "bpstack_1kmvptx_7nct.pdf" # stacked barplots fn
hm.name = "hm-anno-mean-var_nct_notitle.pdf" # heatmaps fn

# format filtered cg datasets
lfilt = list()
lf = list.files()[grepl(".*data.table.filt$", list.files())]
for(f in lf){lfilt[[gsub("\\..*", "", f)]] = fread(f, sep = " ", header = T)}
save(lfilt, file = lfilt.name)

#---------------------------
# quantile mvps (bin method)
#---------------------------
{
  # params
  {
    qiv = seq(0, 1, 0.01) # quantile filter
    qwhich = c(81, 86, 91, 96, 100) # c(17, 18, 19, 20)
    binv = seq(0, 1, 0.01)[1:100] # binned bval mean
  }
  # iter on ncts
  lmvp = list()
  for(t in 1:length(lfilt)){
    ba = lfilt[[t]]
    lci = list(c(), c(), c(), c(), c())
    # iterate on betaval bins
    for(b in binv){
      bf = ba[ba$mean >= b & ba$mean < b + 0.01, ] # filtered on mean bin
      qf = quantile(bf$var, qiv)[qwhich]
      # iterate on ci values
      for(q in 1:length(qf)){
        lci[[q]] = c(lci[[q]], bf$cgid[bf$var > qf[q]])
      }
      message(t, ":", b)
    }
    names(lci) = paste0("ci:", names(qf))
    lmvp[[names(lfilt)[t]]] = lci
  }
  lmvp.bin = lmvp
}
save(lmvp.bin, file = lmvp.bin.name)

#----------------------------------
# quantile mvps (abs cutoff method)
#----------------------------------
lmvp = list()
qiv = seq(0, 1, 0.01) # quantile filter
qwhich = c(81, 86, 91, 96, 100) # c(17, 18, 19, 20)

for(t in 1:length(lfilt)){
  ba = lfilt[[t]]
  lci = list(c(), c(), c(), c(), c())
  qf = quantile(ba$var, qiv)[qwhich]
  for(q in 1:length(qf)){
    lci[[q]] = ba$cgid[ba$var > qf[q]]
  }
  names(lci) = paste0("ci:", names(qf))
  lmvp[[names(lfilt)[t]]] = lci
}

lmvp.abs = lmvp
save(lmvp.abs, file = lmvp.abs.name)

#-------------------------------------------------
# Hypervariable probe types, tissue-specific, etc.
#-------------------------------------------------
# grabs 1k unique tissue-specific hypervar cgs from 2 methods (2k mvps total)
# collapse the lmvp lists
lmvp = list()
# names(lmvp.abs) == names(lmvp.bin)
for(t in names(lmvp.abs)){
  lii = list(); li1 = lmvp.abs[[t]]; li2 = lmvp.bin[[t]]
  for(i in names(li1)){
    lii[[i]] = unique(c(li1[[i]], li2[[i]]))
  }
  lmvp[[t]] = lii
}
save(lmvp, file = lmvp.name)

datlf = lltx = lmvp
txnames = names(datlf)
{
  txmvpl = list()
  fnstem = "filt"
  bpdf1 = bpdf5 = bpdf10 = matrix(nrow = 0, ncol = 3)
  ltx = c(); ltx2 = c()
  mvpcutl = list()
  ptmax = 7 # num shared tx for pan-tissue cgid
  #for(wi in c(6, 5, 4)){
  for(wi in c(5, 4, 3)){
    which.lltx = wi; txmvp = matrix(nrow = 0, ncol = 2)
    for(tx in txnames){
      mvp.id = lltx[[tx]][[which.lltx]]
      txmvp = rbind(txmvp, data.frame(mvp.id, tx, stringsAsFactors = F))
    }
    
    # get pan-tissue probes
    dff = as.data.frame(table(txmvp$mvp.id));
    ptid = dff[dff[,2] == ptmax, 1]; num.ptid = length(ptid)
    if(wi == 5){
      save(ptid, file = "cgid-hypervar99_pantx_7nct.rda")
    }
    if(wi == 4){
      save(ptid, file = "cgid-hypervar95_pantx_7nct.rda")
    }
    if(wi == 3){
      save(ptid, file = "cgid-hypervar90_pantx_7nct.rda")
    }
    
    # filt for tx-specific probes only
    nontxid = c(txmvp[,1][duplicated(txmvp[,1])]) # cgids in >1 tx
    txid = txmvp[,1][!txmvp[,1] %in% nontxid]
    txmvp = txmvp[txmvp[,1] %in% txid,]
    # bp data
    for(tx in txnames){
      mvpidt = lltx[[tx]][[which.lltx]]
      mvpcutl[[tx]] = mvpidt[mvpidt %in% txid]
      datt = table(mvpidt %in% txid)
      nm = matrix(c((as.numeric(datt)[1] - num.ptid), "Non-specific, other", tx,
                    as.numeric(datt)[2], "Tissue-specific", tx,
                    as.numeric(num.ptid), "Pan-tissue", tx), 
                  nrow = 3, byrow = T)
      if(which.lltx == 5){
        bpdf1 = rbind(bpdf1, nm)
      }
      if(which.lltx == 4){
        bpdf5 = rbind(bpdf5, nm)
      }
      if(which.lltx == 3){
        bpdf10 = rbind(bpdf10, nm)
      }
    }
    txmvp = as.data.frame(txmvp, stringsAsFactors = F)
    colnames(txmvp) = c("cgmvp", "txname")
    txmvpl[[length(txmvpl) + 1]] = txmvp
  }
  names(txmvpl) = c("top1.above99.var", "top5.above95.var", "top10.above90.var")
  save(txmvpl, file = paste(txmvpl.name, fnstem, ".rda", sep = ""))
  
}

# args
txmvpl = list()
fnstem = "filt"
bpdf1 = bpdf5 = bpdf10 = matrix(nrow = 0, ncol = 3)
ltx = c(); ltx2 = c()
mvpcutl = list()
ptmax = 7 # num shared tx for pan-tissue cgid
wi = 1
txmvp = matrix(nrow = 0, ncol = 2)
# get all mvp cgids
for(tx in txnames){
  mvp.id = lltx[[tx]][[wi]]
  txmvp = rbind(txmvp, data.frame(mvp.id, tx, stringsAsFactors = F))
}

# get the tx specific cgids
dff = as.data.frame(table(txmvp$mvp.id)) # mvp cgid freq across tissues
cgid.tx = dff$Var1[dff$Freq == 1]  # tx specific

# tx top 1k var cgids, 2 methods
ncgmax = 1000 # max probes by methods
txfiltl = list(); mt = txmvpl[[2]]
for(tx in txnames){
  ltx = lfilt[[tx]]
  # abs method filt
  txmvp.all = lmvp.abs[[tx]][[1]]
  tx.mvpid = txmvp.all[txmvp.all %in% cgid.tx]
  ltxf = ltx[ltx$cgid %in% tx.mvpid,]
  abs.mvp.cgid = ltxf[rev(order(ltxf$var)),]$cgid[1:ncgmax]
  # bin method filt, excluding abs method cgids
  txmvp.all = lmvp.bin[[tx]][[1]]
  tx.mvpid = txmvp.all[txmvp.all %in% cgid.tx & !txmvp.all %in% abs.mvp.cgid]
  ltxf = ltx[ltx$cgid %in% tx.mvpid,]
  bin.mvp.cgid = ltxf[rev(order(ltxf$var)),]$cgid[1:ncgmax]
  # write cgids
  txfiltl[[tx]] = ltx[ltx$cgid %in% c(bin.mvp.cgid, abs.mvp.cgid),]
  message(tx)
}

save(txfiltl, file = txfiltl.name)

#----------------
# make supp table
#----------------

{
  t = "adipose"
  st = as.data.frame(txfiltl[[t]], stringsAsFactors = F)
  for(c in 2:ncol(st)){
    st[,c] = round(as.numeric(st[,c]), 3)
  }
  st = cbind(data.frame(tissue = rep(t, nrow(st)), stringsAsFactors = F),
             st)
  cgf = cga[cga$Name %in% st$cgid,]
  cgf = cgf[order(match(cgf$Name, st$cgid)),]
  identical(cgf$Name, st$cgid)
  sta = cbind(st, cgf[,c(1, 2, 3, 9, 14, 15, 20, 21, 22)])
  
  for(t in names(txfiltl)[2:length(txfiltl)]){
    st = as.data.frame(txfiltl[[t]], stringsAsFactors = F)
    for(c in 2:ncol(st)){
      st[,c] = round(as.numeric(st[,c]), 3)
    }
    st = cbind(data.frame(tissue = rep(t, nrow(st)), stringsAsFactors = F),
               st)
    cgf = cga[cga$Name %in% st$cgid,]
    cgf = cgf[order(match(cgf$Name, st$cgid)),]
    if(identical(cgf$Name, st$cgid)){
      st = cbind(st, cgf[,c(1, 2, 3, 9, 14, 15, 20, 21, 22)])
    } else{
      message("sets not matched!")
    }
    sta = rbind(sta, st)
  }
}

save(sta, file = paste0(sta.name, ".rda"))
write.csv(sta, file = paste0(sta.name, ".csv"))


#------------------------------------------------
# get st pan-tissue hypovar cg (quantile method)
#------------------------------------------------
# get bottom 10th quantile per tissue
i = 1
di = lfilt[[i]]
which.cg = di$var < quantile(di$var, seq(0, 1, 0.1))[2]
pan.cghypo = di[which.cg,]$cgid # begin pan-tissue cg id list

for(i in 2:length(lfilt)){
  di = lfilt[[i]]
  which.cg = di$var < quantile(di$var, seq(0, 1, 0.1))[2]
  pan.cghypo = intersect(pan.cghypo, di[which.cg,]$cgid)
}
st = cga[cga$Name %in% pan.cghypo, c(1, 2, 3, 9, 14, 15, 20, 22)]

write.csv(st, file = paste0(st.name, ".csv"))
save(st, file = paste0(st.name, ".rda"))

#-------------
# violin plots
#-------------
txfiltl <- get(load(txfiltl.name))
txnames = names(txfiltl)
{
  bpdf.mean = matrix(nrow = 0, ncol = 3)
  bpdf.var = matrix(nrow = 0, ncol = 3)
  # tx colors
  tl = c("blood", "buccal", "brain", "sperm", "nasal", "adipose", "liver" ) # index
  tc = c("red", "orange", "purple", "blue", "green", "brown", "forestgreen") # color
  for(t in 1:length(txnames)){
    tx = txnames[t]; datt = txfiltl[[tx]]
    bpdf.mean = rbind(bpdf.mean, data.frame(datt$mean, tc[tl == tx], tx, stringsAsFactors = F))
    bpdf.var = rbind(bpdf.var, data.frame(datt$var, tc[tl == tx], tx, stringsAsFactors = F))
  }
  bpdf.mean = as.data.frame(bpdf.mean, stringsAsFactors = F)
  bpdf.var = as.data.frame(bpdf.var, stringsAsFactors = F)
  bpdf.mean[,1] = as.numeric(as.character(bpdf.mean[,1]))
  bpdf.var[,1] = as.numeric(as.character(bpdf.var[,1]))
  colnames(bpdf.mean) = c("mean", "col", "tissue")
  colnames(bpdf.var) = c("var", "col", "tissue")
  
  bpdf.mean$tissue = factor(bpdf.mean$tissue, levels = unique(bpdf.mean$tissue))
  ordert = order(match(tl, levels(bpdf.mean$tissue)))

  # version 1 (composite plot)
  {
    p1 = ggplot(bpdf.mean, aes(x = tissue, y = mean, fill = tissue)) + 
      geom_violin(trim = F, show.legend = F) +
      scale_fill_manual(values = tc[ordert]) + theme_bw() + 
      theme(axis.text.x = element_text(angle = 90)) + 
      ggtitle("") + xlab("") + ylab("Mean")
    
    bpdf.var$tissue = factor(bpdf.var$tissue, levels = unique(bpdf.var$tissue))
    ordert = order(match(tl, levels(bpdf.var$tissue)))
    
    p2 = ggplot(bpdf.var, aes(x = tissue, y = var, fill = tissue)) + 
      geom_violin(trim = F, show.legend = F) +
      scale_fill_manual(values = tc[ordert]) + theme_bw() + 
      theme(axis.text.x = element_text(angle = 90)) + 
      ggtitle("") + xlab("") + ylab("Variance")
    
    pdf(vp.name, 5, 5)
    grid.arrange(p1, p2, ncol = 1, 
                 bottom = "Tissue")
    dev.off()
  }
  
  # version 2 (separate plots)
  {
    p1 = ggplot(bpdf.mean, aes(x = tissue, y = mean, fill = tissue)) + 
      geom_violin(trim = F, show.legend = F) +
      scale_fill_manual(values = tc[ordert]) + theme_bw() + 
      theme(axis.text.x = element_text(angle = 90)) +
      xlab("Tissue") + ylab("Mean") +
      scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, 0.5))
    
    pdf("vp-mean_hivar_7nct.pdf", 3, 2)
    p1
    dev.off()
    
    bpdf.var$tissue = factor(bpdf.var$tissue, levels = unique(bpdf.var$tissue))
    ordert = order(match(tl, levels(bpdf.var$tissue)))
    
    p2 = ggplot(bpdf.var, aes(x = tissue, y = var, fill = tissue)) + 
      geom_violin(trim = F, show.legend = F) +
      scale_fill_manual(values = tc[ordert]) + theme_bw() + 
      theme(axis.text.x = element_text(angle = 90)) + 
      xlab("Tissue") + ylab("Variance") +
      scale_y_continuous(limits = c(0, 0.12), breaks = seq(0, 0.12, 0.06))
    
    pdf("vp-var_hivar_7nct.pdf", 3, 2)
    p2
    dev.off()
    
    
  }
  
}

pdf(vp.name, 5, 5)
grid.arrange(p1, ncol = 1, 
             bottom = "Tissue", left = "  Mean DNAm (Beta-value)")
dev.off()

#---------------------
# add modified cg anno
#---------------------
prom.stat = grepl("TSS|5'", cga$UCSC_RefGene_Group)
body.stat = grepl("Body|Exon|3'", cga$UCSC_RefGene_Group)
cga$gene.type = ifelse(cga$UCSC_RefGene_Name=="", "intergenic",
                       ifelse(prom.stat & body.stat, "intragenic_promoter-body",
                              ifelse(prom.stat & !body.stat, "intragenic_promoter", "intragenic_body")))
cga$isl.type = ifelse(cga$Relation_to_Island=="OpenSea", "interisland_opensea", 
                      ifelse(cga$Relation_to_Island=="Island", "intraisland_main", "intraisland_other"))
table(cga$isl.type, cga$gene.type)
cga$type.composite = paste0(cga$isl.type,";",cga$gene.type)

cga$annocat = ifelse(cga$Relation_to_Island=="OpenSea" & cga$UCSC_RefGene_Name=="", "interisland_intergenic",
                     ifelse(!(cga$Relation_to_Island=="OpenSea") & cga$UCSC_RefGene_Name=="", "intraisland_intergenic",
                            ifelse(cga$Relation_to_Island=="OpenSea" & !(cga$UCSC_RefGene_Name==""), "interisland_intragenic",
                                   ifelse(!(cga$Relation_to_Island=="OpenSea" | cga$UCSC_RefGene_Name==""), "intraisland_intragenic", "NA"))))

prom.stat = grepl("TSS|5'", cga$UCSC_RefGene_Group)
body.stat = grepl("Body|Exon|3'", cga$UCSC_RefGene_Group)
cga$gene.relation = ifelse(cga$UCSC_RefGene_Name=="", "intergenic",
                           ifelse(prom.stat & body.stat, "promoter;body",
                                  ifelse(prom.stat & !body.stat, "promoter", "body")))

#------------------------
# barplots, mapping freqs
#------------------------
txnames = names(txfiltl)
mainplotwidth = 8; mainplotheight = 8
legendpos <- "top"
# anno plots: isl/gene relation
{
  # top 1k mvp, tissue-specific
  bpdf = matrix(nrow = 0, ncol = 3)
  which.lltx = 6
  for(ti in 1:length(txnames)){
    txname = txnames[ti]; mvpidt = txfiltl[[txname]]$cgid
    tdati = as.data.frame(table(cga[rownames(cga) %in% mvpidt,]$annocat))
    tdati$txname = txname; bpdf = rbind(bpdf, tdati)
  }
  colnames(bpdf) = c("Genomic Region", "Count", "Tissue")
  bpdf$Tissue = factor(bpdf$Tissue, levels=txnames)
  
  # get compact label from region var
  bpdf$label <- "NA"
  bpdf[bpdf$`Genomic Region` == "interisland_intergenic",]$label <- "Osea;Nogene"
  bpdf[bpdf$`Genomic Region` == "interisland_intragenic",]$label <- "Osea;Gene"
  bpdf[bpdf$`Genomic Region` == "intraisland_intergenic",]$label <- "Isld;Nogene"
  bpdf[bpdf$`Genomic Region` == "intraisland_intragenic",]$label <- "Isld;Gene"
  
  # new plot, compact legend format
  p1 = ggplot(bpdf, aes(x = Tissue, y = Count, fill = label)) +
    geom_bar(stat = "identity") + theme_bw() +
    theme(axis.text.x = element_text(angle = 90),
          legend.title = element_blank(),
          legend.position = legendpos,
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
}
# anno plots: gene relation
{
  bpdf = matrix(nrow = 0, ncol = 3)
  which.lltx = 6
  for(ti in 1:length(txnames)){
    txname = txnames[ti]; mvpidt = txfiltl[[txname]]$cgid
    tdati = as.data.frame(table(cga[rownames(cga) %in% mvpidt,]$gene.relation))
    tdati$txname = txname; bpdf = rbind(bpdf, tdati)
  }
  colnames(bpdf) = c("Gene Region", "Count", "Tissue")
  bpdf$Tissue = factor(bpdf$Tissue, levels=txnames)
  
  # new compact label from region var
  bpdf$label <- "NA"
  bpdf[bpdf$`Gene Region` == "body",]$label <- "Body"
  bpdf[bpdf$`Gene Region` == "promoter",]$label <- "Prom"
  bpdf[bpdf$`Gene Region` == "intergenic",]$label <- "Nogene"
  bpdf[bpdf$`Gene Region` == "promoter;body",]$label <- "Prom;Body"
  
  # new plot with compact legend
  p2 = ggplot(bpdf, aes(x = Tissue, y = Count, fill = label)) +
    geom_bar(stat = "identity") + theme_bw() +
    theme(axis.text.x = element_text(angle = 90),
          legend.position = legendpos,
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) + 
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
}
# anno plots: isl relation
{
  # top 1k mvp, tissue-specific
  bpdf = matrix(nrow = 0, ncol = 3)
  which.lltx = 6
  for(ti in 1:length(txnames)){
    txname = txnames[ti]; mvpidt = txfiltl[[txname]]$cgid
    tdati = as.data.frame(table(cga[rownames(cga) %in% mvpidt,]$Relation_to_Island))
    tdati$txname = txname; bpdf = rbind(bpdf, tdati)
  }
  colnames(bpdf) = c("Island Region", "Count", "Tissue")
  bpdf$Tissue = factor(bpdf$Tissue, levels=txnames)
  
  # new compact label from region var
  bpdf$label <- "NA"
  bpdf[bpdf$`Island Region` == "OpenSea",]$label <- "OSea"
  bpdf[bpdf$`Island Region` == "Island",]$label <- "Isld"
  bpdf[bpdf$`Island Region` == "N_Shelf",]$label <- "Nshl"
  bpdf[bpdf$`Island Region` == "N_Shore",]$label <- "Nshr"
  bpdf[bpdf$`Island Region` == "S_Shelf",]$label <- "Sshl"
  bpdf[bpdf$`Island Region` == "S_Shore",]$label <- "Sshr"
  
  # new plot with compact legend
  p3 = ggplot(bpdf, aes(x = Tissue, y = Count, fill = label)) +
    geom_bar(stat = "identity") + theme_bw() +
    theme(axis.text.x = element_text(angle = 90),
          legend.position = legendpos,
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) + 
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
}
# version 1 (vertical arrangement)
pdf("bp-region_hi_nct7.pdf", 6, 8)
grid.arrange(p1, p2, p3, ncol = 1,
             bottom="Tissue", left="Number of\nProbes")
dev.off()
# version 2 (horizontal arrangement)
pdf("bp-region_hi_nct7.pdf", 7.2, 2.5)
grid.arrange(p1, p2, p3, ncol = 3,
             bottom = "Tissue", 
             left = "Number of\nProbes")
dev.off()

#---------------------------------
# hm meth summaries by region, nct
#---------------------------------
mincg = 2
txnames = names(txfiltl)
annom = matrix(nrow = 0, ncol = length(txnames))
utype = unique(cga$type.composite)
hmma.mean = hmma.var = hmma.size = matrix(nrow = 0, ncol = length(txnames))

for(a in utype){
  cgidv = rownames(cga[cga$type.composite == a,])
  rdat.mean = rdat.var = rdat.size = c()
  for(tx in txnames){
    datx = txfiltl[[tx]]; mvpt = datx$cgid
    cgint = intersect(mvpt, cgidv); cgint.size = length(cgint)
    if(length(cgint) >= mincg){
      rdat.mean = c(rdat.mean, mean(datx[datx$cgid %in% cgint,]$mean))
      rdat.var = c(rdat.var, var(datx[datx$cgid %in% cgint,]$mean))
      rdat.size = c(rdat.size, cgint.size)
    } else{
      rdat.mean = c(rdat.mean, "NA"); rdat.var = c(rdat.var, "NA")
    }
    message(tx)
  }
  hmma.mean = rbind(hmma.mean, rdat.mean)
  hmma.var = rbind(hmma.var, rdat.var)
  hmma.size = rbind(hmma.size, rdat.size)
  message(a)
}

rownames(hmma.mean) = rownames(hmma.var) = rownames(hmma.size) = utype
colnames(hmma.mean) = colnames(hmma.var) = colnames(hmma.size) = txnames
class(hmma.mean) = class(hmma.var) = class(hmma.size) = "numeric"

hdx = hdv = matrix(nrow = 0, ncol = 4)
for(r in rownames(hmma.mean)){
  for(c in colnames(hmma.mean)){
    hdx = rbind(hdx, matrix(c(hmma.mean[r, c], c, r, hmma.size[r, c]), nrow = 1))
    hdv = rbind(hdv, matrix(c(hmma.var[r, c], c, r, hmma.size[r, c]), nrow = 1))
  }
}

hdx = as.data.frame(hdx, stringsAsFactors = F)
hdv = as.data.frame(hdv, stringsAsFactors = F)
colnames(hdx) = c("mean", "tissue", "anno", "size")
colnames(hdv) = c("var", "tissue", "anno", "size")
hdx$mean = as.numeric(hdx$mean); hdx$size = as.numeric(hdx$size)
hdv$var = as.numeric(hdv$var); hdv$size = as.numeric(hdv$size)

hm.mean = ggplot(hdx, aes(tissue, anno)) +
  geom_tile(aes(fill = mean)) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5,
                       limits = c(0, 1)) +
  geom_text(aes(label = size), color = "black") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Mean") + ylab("") + xlab("")

hm.var = ggplot(hdv, aes(tissue, anno)) +
  geom_tile(aes(fill = var)) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0.04, limits = c(0, 0.142)) +
  geom_text(aes(label = size), color = "black") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  xlab("") + ylab("") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ggtitle("Variance")

pdf(hm.name, 12, 5)
grid.arrange(hm.mean, hm.var, 
             layout_matrix = matrix(c(rep(1, 7), rep(2, 4)), nrow = 1),
             bottom = "Tissue", left = "Annotation/Region Type")
dev.off()