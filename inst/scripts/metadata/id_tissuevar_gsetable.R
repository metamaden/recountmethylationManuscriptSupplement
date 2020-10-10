
# filter on sample type
tissue.terms <- c("sperm","liver","kidney","breast","muscle",
                  "prostate","skin","bone","stomach","esophagus")

tdf <- data.frame(tissue="",gse="")

for(t in 1:length(tissue.terms)){
  gt <- unique(xdf[grepl(tissue.terms[t],xdf$sample_type),]$gse)
  for(g in 1:length(gt)){
    tdf <- rbind(tdf, data.frame(tissue=tissue.terms[t],
                                 gse=gt[g]))
  }
  message(t)
}

write.csv(tdf, file="gse_tissuevar_gsetable.csv")

# filter on additional terms
xdf[xdf$anatomic_location=="skin",]