#!/usr/bin/env R

# Summarize probes with high tissue-specific autosomal 
# DNAm variance across 7 noncancer tissues (Table S8).

#----------
# load data
#----------
pkgname <- "recountmethylationManuscriptSupplement"
nct7.dir <- system.file("nct7", "pcadata", package = pkgname) 
fn <- "lfilt_highvar-nc7.rda"
txfiltl <- get(load(file.path(nct7.dir, fn)))

#---------------
# get table data
#---------------
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

#save(sta, file = paste0(sta.name, ".rda"))
#write.csv(sta, file = paste0(sta.name, ".csv"))