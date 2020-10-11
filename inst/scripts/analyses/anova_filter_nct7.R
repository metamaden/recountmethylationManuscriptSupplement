#!/usr/bin/env R

# Author: Sean Maden
# Perform ANOVAs on batch corrected noob-normalized autosomal DNAm in 
# 7 noncancer tissues. Note the ANOVA results tables are available from
# the data server at https://recount.bio/data/recountmethylation_manuscript_supplement/data/

library(data.table)

#----------
# load data
#----------
# get minfi-formatted probe annotations
cga <- getAnnotation(get(data("MsetEx")))

# notes: no pval for residuals, and no predsex for sperm
lbc <- list(); adn <- "anova_results"
txl <- c("adipose", "blood", "brain", "buccal", "liver", "nasal", "sperm")

# read anova results
for(t in txl){
  prebct <- fread(paste(adn, paste0("prebc.anova.table.", t), sep = "/"), sep = " ", data.table = F)
  postbct <- fread(paste(adn, paste0("postbc.anova.table.", t), sep = "/"), data.table = F)
  if(!t == "sperm"){
    {
      # fix colnames
      pref <- prebct[,c(1:(ncol(prebct)-1))]
      colnames(pref) <- c(colnames(prebct)[1:9], colnames(prebct)[11:ncol(prebct)])
      postf <- postbct[,c(1:(ncol(postbct)-1))]
      colnames(postf) <- c(colnames(postbct)[1:9], colnames(postbct)[11:ncol(postbct)])
      # get padj vars
      ncn <- c("predage.padj", "predsex.padj", "predcd8t.padj", "predcd4t.padj", "prednk.padj", 
              "predbcell.padj", "predmono.padj", "predgran.padj")
      for(c in 2:9){
        pref$newcol <- p.adjust(pref[,c]); colnames(pref)[ncol(pref)] <- ncn[c-1]
        postf$newcol <- p.adjust(postf[,c]); colnames(postf)[ncol(postf)] <- ncn[c-1]
      }
    }
  } else{
    # fix colnames
    which.cols.pre <- c(1:2, 4:9, 11, 13:20, 22:28)
    which.cols.post <- c(1:8, 10:25)
    pref <- prebct[, c(1:24)]; colnames(pref) <- colnames(prebct)[which.cols.pre]
    postf <- postbct[, c(1:24)]; colnames(postf) <- colnames(postbct)[which.cols.post]
    # get padj vars
    ncn <- c("predage.padj", "predcd8t.padj", "predcd4t.padj", "prednk.padj", 
            "predbcell.padj", "predmono.padj", "predgran.padj")
    for(c in 2:8){
      pref$newcol <- p.adjust(pref[,c]); colnames(pref)[ncol(pref)] <- ncn[c-1]
      postf$newcol <- p.adjust(postf[,c]); colnames(postf)[ncol(postf)] <- ncn[c-1]
    }
  }
  lbc[[t]] <- list("prebc" = pref, "postbc" = postf)
  message(t)
}
save(lbc, file = "ldat-anovabc_nct.rda")

#-----------------------
# read cg summary tables
#-----------------------
dn <- "cgall"
lf <- list.files("cgall")
cg.adipose <- fread(paste(dn, lf[grepl("adipose", lf)], sep = "/"), sep = " ", header = T, data.table = F)
cg.blood <- fread(paste(dn, lf[grepl("blood", lf)], sep = "/"), sep = " ", header = T, data.table = F)
cg.brain <- fread(paste(dn, lf[grepl("brain", lf)], sep = "/"), sep = " ", header = T, data.table = F)
cg.liver <- fread(paste(dn, lf[grepl("liver", lf)], sep = "/"), sep = " ", header = T, data.table = F)
cg.buccal <- fread(paste(dn, lf[grepl("buccal", lf)], sep = "/"), sep = " ", header = T, data.table = F)
cg.nasal <- fread(paste(dn, lf[grepl("nasal", lf)], sep = "/"), sep = " ", header = T, data.table = F)
cg.sperm <- fread(paste(dn, lf[grepl("sperm", lf)], sep = "/"), sep = " ", header = T, data.table = F)
datl <- list(cg.blood, cg.brain, cg.adipose, cg.liver, cg.nasal, cg.buccal, cg.sperm); 
names(datl) <- c("blood", "brain", "adipose", "liver", "nasal", "buccal", "sperm")

#-------------------------
# apply anova and chr filt
#-------------------------
# filter params for anova vars
pfilt <- 1e-3 # pvalue of variance
vpercfilt <- 10 # percent variance
cgid.chrfilt <- cga[cga$chr %in% c("chrX", "chrY"),]$Name

# filter on anova results
bpdf <- matrix(nrow = 0, ncol = 3)
for(t in names(lbc)){
  postbct <- lbc[[t]]$postbc
  cgdatt <- datl[[t]]
  # chr filt
  cgdatt <- cgdatt[!cgdatt$cgid %in% cgid.chrfilt,]
  # anova filt
  filt <- (postbct$predage.padj < pfilt & postbct$predagevarperc >= vpercfilt) |
    (postbct$predcd8t.padj < pfilt & postbct$cd8tvarperc >= vpercfilt) | 
    (postbct$predcd4t.padj < pfilt & postbct$cd4tvarperc >= vpercfilt) |
    (postbct$prednk.padj < pfilt & postbct$nkvarperc >= vpercfilt) | 
    (postbct$predbcell.padj < pfilt & postbct$bcellvarperc >= vpercfilt) | 
    (postbct$predmono.padj < pfilt & postbct$monovarperc >= vpercfilt) |
    (postbct$predgran.padj < pfilt & postbct$granvarperc >= vpercfilt)
  if(!t == "sperm"){
    filt <- filt | (postbct$predsex.padj < pfilt & postbct$predsexvarperc >= vpercfilt)
  }
  # apply filt
  rmcg <- postbct[filt,]
  cgdat.filt <- cgdatt[!cgdatt$cgid %in% rmcg$cgid,]
  # write new data
  tf <- table(filt)
  nm <- matrix(c(tf[2], "REMOVED", t, tf[1], "RETAINED", t), ncol = 3, byrow=T)
  bpdf <- rbind(bpdf, nm)
  fwrite(cgdat.filt, paste(t, "data.table.filt", sep = "."), sep = " ")
}
colnames(bpdf) <- c("nprobes", "type", "tissue")
bpdf <- as.data.frame(bpdf, stringsAsFactors = F)
bpdf$nprobes <- as.numeric(bpdf$nprobes)
tv <- unique(bpdf$tissue); bf <- bpdf[bpdf$type =="REMOVED",]
bpdf$tissue <- factor(bpdf$tissue, levels = bf[order(bf$nprobes),]$tissue)

# save bpdf for anova filters plot
fwrite(bpdf, "anova-filt_nct.table", sep = " ")

#----------------------------------
# format data tables as list object
#----------------------------------
lfilt <- list()
lf <- list.files()[grepl(".*data.table.filt$", list.files())]
for(f in lf){lfilt[[gsub("\\..*", "", f)]] <- fread(f, sep = " ", header = T)}
save(lfilt, file = lfilt.name)
