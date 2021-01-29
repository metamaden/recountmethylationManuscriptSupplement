#!/usr/bin/env R

# Author: Sean Maden
# Summarize probe and genome region mappings at probes with low variances 
# and mean ranges across 7 noncancer tissues (Table S7).

library(minfiData)

#----------
# load data
#----------
# get minfi-formatted probe annotations
cga <- getAnnotation(get(data("MsetEx")))
# get probe ids
pkgname <- "recountmethylationManuscriptSupplement"
nct7.dir <- system.file("extdata", "nct7", package = pkgname)
id.fn <- "cgids_lowvar-lowmean_nct7.rda"
ptid <- get(load(file.path(nct7.dir, id.fn)))
# get probe summary statistics
lfilt <- get(load(file.path(nct7.dir, "lfilt_cgtables_nct.rda")))

#---------------
# make new table
#---------------
# get anno and match data
cnames <- c("Name", "chr", "pos", "strand", "Type",	
            "Islands_Name", "Relation_to_Island",	
            "UCSC_RefGene_Name", "UCSC_RefGene_Accession",	
            "UCSC_RefGene_Group")
which.cn <- which(colnames(cga) %in% cnames)
cgf <- cga[cga$Name %in% ptid, which.cn]

#--------------------------------
# append dnam variances by tissue
#--------------------------------
mt <- data.frame(cgid = cgf$Name, stringsAsFactors = FALSE)
for(t in names(lfilt)){
  lf <- lfilt[[t]];df <- lf[lf$cgid %in% cgf$Name,]
  df <- df[order(match(df$cgid, cgf$Name)),]
  cond <- identical(df$cgid, cgf$Name)
  if(cond){
    dft <- df[,c("cgid", "var", "mean")]
    colnames(dft) <- c("cgid", paste0(t,c(".variance", ".mean")))
    # significant figures
    dft[,2] <- format(dft[,2], scientific = TRUE, digits = 3)
    dft[,3] <- format(dft[,3], scientific = TRUE, digits = 3)
    mt <- cbind(mt, dft[,c(2:3)])
  }else{stop("Error matching cgids for tissue ", t)};message(t)
}

mt$max.interval.means <- apply(mt, 1, function(x){
  mtf <- x[grepl("mean", colnames(mt))]
  mintv <- unlist(lapply(unique(gsub("\\..*", "", names(mtf))), function(x){
    gfx <- grepl(x, names(mtf));m1 <- mtf[gfx]; m2v <- mtf[!gfx]
    diffvx <- rep(as.numeric(m1), length(m2v)) - as.numeric(m2v)
    return(abs(diffvx))}))
  return(max(mintv))})

cond <- identical(mt$cgid, cgf$Name)
if(cond){cgf <- cbind(cgf, mt[,c(2:ncol(mt))])}else{
  stop("Error matching cgid's for mt, cgf.")}

#------
# save
#------
st <- cgf; st.fn <- "table-s7_recurr-lowvar-cg-anno_7nct"

# no caption
save(st, file = paste0(st.fn, ".rda"))
write.csv(st, file = paste0(st.fn, ".csv"), row.names = FALSE)

# with caption
cap.txt <- paste0("Table S7. Characteristics of DNAm array ",
                    "CpG probes (rows) with recurrent low Beta-value ",
                    "variances and mean intervals across 7 non-cancer ",
                  "tissues. Column 1 is the probe ID. Column 2 is the ",
                  "chromosome. Column 3 ",
                  "is genome coordinate. Column 4 is the DNA strand. ",
                  "Column 5 the assay type. Column 6 is the CpG island name. ",
                  "Column 7 is the CpG island region type. Column 8 is the ",
                  "gene ID(s). Column 9 is the gene accession(s). Column 10 ",
                  "is the gene region group(s). Columns 11-24 are the ",
                  "Beta-value variances and means across samples by tissue ",
                  "(tissues: adipose; blood; brain; buccal; liver; nasal; ",
                  "and sperm). Column 25 is the maximum absolute mean ",
                  "interval across all pairwise comparisons of tissue means.")
file.conn <- file(paste0(st.fn, ".csv")); writeLines(cap.txt, file.conn)
write.table(st, file = paste0(st.fn, ".csv"), row.names = FALSE, 
            append = TRUE, sep = ",")
(tissues: adipose, blood, brain, buccal, liver, nasal, and sperm)