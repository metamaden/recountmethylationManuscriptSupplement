#!/usr/bin/env R

# Make acronyms key, Table 1.

#-----------------------
# define table cell data
#-----------------------
t1 <- data.frame(Term = c("GEO", "GSE", "GSM", "CpG", "DNAm", "HM450K"),
                 Meaning = c("Gene expression omnibus", 
                             "Study accession number",
                             "Sample accession number",
                             "Cytosine-guanine dinucleotide",
                             "DNA methylation","HumanMethylome 450K"),
                 Description = c(paste0("The largest public database of ",
                                        "published array data."),
                                paste0("Unique identifier for a study ",
                                       "record that includes a platform, ",
                                       "set of sample records, and ",
                                       "supplemental matrices containing ",
                                       "assay data."),
                                paste0("Unique identifier for a sample record ",
                                       "that includes sample-specific metadata ",
                                       "and may include supplemental sample ",
                                       "datasets."),
                                paste0("Dinucleotide sequence, or locus, ",
                                       "consisting of a cytosine followed by ",
                                       "a guanine."),
                                paste0("The presence of a nucleotide-bound ",
                                       "methyl group, typically at the 5' ",
                                       "cytosine position in a CpG locus."),
                                paste0("Popular array platform, manufactured ",
                                       "by Illumina, that uses BeadArray ",
                                       "technology to probe DNA methylation ",
                                       "at roughly 480,000 CpG loci.")),
                 stringsAsFactors = FALSE)

# caption text 
captxt <- "Meanings and definitions of frequently used terms and acronyms."

#-----------------
# save/print table
#-----------------
# save table with caption
t1.fn <- "table1.csv"; cat(paste0("Table 1. ", captxt, "\n"), file = t1.fn)
write.table(t1, file = t1.fn, sep = ",", row.names = FALSE, append = TRUE)

# table latex script
ld <- list();
ld[[length(ld) + 1]] <- "\\begin{table}[h]"
ld[[length(ld) + 1]] <- "\\begin{center}"
ld[[length(ld) + 1]] <- paste0("\\caption{",captxt,"}")
ld[[length(ld) + 1]] <- "\\begin{tabular}{ | c | p{4cm}| p{9cm} | } "
ld[[length(ld) + 1]]<-paste0("\\hline")
cv <- colnames(t1);rval <- paste(cv[1], cv[2], cv[3], sep = " & ")
ld[[length(ld) + 1]] <- paste0(rval, "\\\\"); ld[[length(ld) + 1]] <- "\\hline"
for(r in seq(nrow(t1))){
  rval <- paste(t1[r,1], t1[r,2], t1[r,3], sep = " & ")
  ld[[length(ld) + 1]] <- paste0(rval, "\\\\")
  ld[[length(ld) + 1]] <- "\\hline"}
ld[[length(ld) + 1]] <- "\\end{tabular}"
ld[[length(ld) + 1]] <- "\\end{center}"
ld[[length(ld) + 1]] <- "\\end{table}"

# print latex script lines
for(i in seq(length(ld))){print(ld[[i]])}

# save latex script lines
file.conn <- file("table1_latex.txt")
writeLines(unlist(ld), file.conn)
close(file.conn)




