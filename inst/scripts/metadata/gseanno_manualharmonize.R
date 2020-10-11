#!/usr/bin/env R

# Author: Sean Maden
# This script describes the coercion of mined metadata into shared term 
# labels for sample type, disease state, gender, age, anatomic location, 
# etc. As input, this script loads the list of study-wise SOFT metadata 
# tables downloaded from GEO (see the script "make_gse_annolist.R" for 
# details). For info about acquiring the SOFT metadata, see the repo
# "recount-methylation-server".

#----------
# load data
#----------
tgse.list <- get(load("datafiles/geo_gse-atables_list.rda"))
gat <- matrix(nrow = 0, ncol = 8)
colnames(gat) <- c("gsm","gse","sample_type","disease_state",
                   "gender","age","anatomic_location","misc")

#---------------------
# process study tables
#---------------------

# GSE126017
i=1
ti <- as.data.frame(tgse.list[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(as.character(ti$gsm[r]),
                             as.character(ti$gse[r]),
                             as.character(ti$`sample type`[r]),
                             as.character(ti$`disease state`[r]),
                             as.character(ti$gender[r]),
                             "NA", "NA", "NA"),nrow=1))
  message(r)
}

# GSE106360
i=2
ti <- as.data.frame(tgse.list[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(as.character(ti$gsm[r]),
                             as.character(ti$gse[r]),
                             as.character(ti$`sample type`[r]),
                             as.character(ti$`disease state`[r]),
                             as.character(ti$gender[r]),
                             "NA", "NA", "NA"),nrow=1))
  message(r)
}

# GSE111428
i=3
ti <- as.data.frame(tgse.list[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(as.character(ti$gsm[r]), # GSM ID
                             as.character(ti$gse[r]), # GSE ID
                             as.character(ti$tumour[r]), # Sample type
                             "cancer", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             as.character(ti$location[r]), # Anatomic Location
                             paste0("cell_line_id:",ti$`cell line id`[r])), 
                           nrow=1))
  message(r)
}

# GSE110780
i=4
ti <- as.data.frame(tgse.list[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(as.character(ti$gsm[r]),
                             as.character(ti$gse[r]),
                             "NA",
                             as.character(ti$dx[r]),
                             as.character(ti$Sex[r]),
                             as.character(ti$age[r]),
                             "NA",
                             paste0("race:",ti$race[r],";",
                                    "treatment:",ti$treatment[r],";",
                                    "group:",ti$group[r],";",
                                    "smoking_status:",ti$`smoking status`[r])),
                           nrow=1))
  message(r)
}

# GSE110778
i=5
ti <- as.data.frame(tgse.list[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(as.character(ti$gsm[r]),
                             as.character(ti$gse[r]),
                             "NA",
                             as.character(ti$dx[r]),
                             as.character(ti$Sex[r]),
                             as.character(ti$age[r]),
                             "NA",
                             paste0("race:",ti$race[r],";",
                                    "ethnicity:",ti$ethnicity[r],";",
                                    "treatment:",ti$treatment[r],";",
                                    "smoking_status:",ti$`smoking status`[r])),
                           nrow=1))
  message(r)
}

# GSE110776
i=6
ti <- as.data.frame(tgse.list[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(as.character(ti$gsm[r]),
                             as.character(ti$gse[r]),
                             "NA",
                             "NA",
                             as.character(ti$Sex[r]),
                             as.character(ti$age[r]),
                             "NA",
                             paste0("race:",ti$race[r],";",
                                    "ethnicity:",ti$ethnicity[r],";",
                                    "group:",ti$group[r],";",
                                    "treatment:",ti$treatment[r],";",
                                    "smoking_status:",ti$`smoking status`[r])),
                           nrow=1))
  message(r)
}

# GSE83944
i=7
ti <- as.data.frame(tgse.list[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "NA", # Disease state
                             ti$Sex[r], # Gender
                             ti$age[r], # Age
                             "NA", # Anatomic Location
                             "NA"), # misc
                           nrow=1))
  message(r)
}

# GSE75704
i=8
ti <- as.data.frame(tgse.list[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tssue[r], # Sample type
                             ti$`subject group`[r], # Disease state
                             ti$gender[r], # Gender
                             ti$`age at death`[r], # Age
                             "NA", # Anatomic Location
                             paste0("age_info:age_at_death;",
                                    "subject_id:", ti$`subject id`[r],";",
                                    "subject_status:", ti$`subject status`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE108423
i=9
ti <- as.data.frame(tgse.list[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "NA", # Sample type
                             "NA", # Disease state
                             ti$Sex[r], # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             "NA"),# misc
                           nrow=1))
  message(r)
}

# study 10: GSE121633
i=10
ti <- as.data.frame(tgse.list[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "NA", # Sample type
                             "NA", # Disease state
                             ti$Sex[r], # Gender
                             ti$age[r], # Age
                             "NA", # Anatomic Location
                             "NA"),# misc
                           nrow=1))
  message(r)
}

# resort by size, and continue
sv <- unlist(lapply(tgse.list,nrow))
tls <- tgse.list[rev(order(sv))]

# GSE51032
i=1
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("material:",ti$material[r])),# misc
                           nrow=1))
  message(r)
}

# GSE90496
i=2
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "NA", # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("material:",ti$material[r])),# misc
                           nrow=1))
  message(r)
}

# GSE105018
i=3
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "NA", # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("twinid:",ti$twinid[r],";",
                                    "familyid:",ti$familyid[r],";",
                                    "zygosity:",ti$zygosity[r]
                             )),# misc
                           nrow=1))
  message(r)
}

# GSE109379
i=4
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "cancer", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("material:",ti$material[r])),# misc
                           nrow=1))
  message(r)
}

# GSE68379
i=5
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`primary histology`[r], # Sample type
                             "cancer", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             ti$`primary site`[r], # Anatomic Location
                             paste0("cell_line:",ti$`cell line`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE51032
i=6
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "NA", # Sample type
                             "NA", # Disease state
                             ti$gender[r], # Gender
                             ti$age[r], # Age
                             "NA", # Anatomic Location
                             paste0("cancer_type:",ti$`cancer type (icd-10)`[r],";",
                                    "time_to_diagnosis:",ti$`time to diagnosis`[r],";",
                                    "age_at_menarche:",ti$`age at menarche`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE85218
i=7
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "cancer", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("subgroup:",ti$subgroup[r],";",
                                    "type:",ti$type[r])),# misc
                           nrow=1))
  message(r)
}

# GSE85212
i=8
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "cancer", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("subgroup:",ti$subgroup[r],";",
                                    "type:",ti$type[r])),# misc
                           nrow=1))
  message(r)
}

# GSE87571
i=9
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             ti$`disease state`[r], # Disease state
                             ti$gender[r], # Gender
                             ti$age[r], # Age
                             "NA", # Anatomic Location
                             "NA"),# misc
                           nrow=1))
  message(r)
}

# GSE74193
i=10
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "NA", # Sample type
                             ti$group[r], # Disease state
                             ti$`sex (clinical gender)`[r], # Gender
                             ti$`age (in years)`[r], # Age
                             "NA", # Anatomic Location
                             paste0("predicted_sex:",ti$`predictedsex (sex predicted by sex chromosome probes)`[r],";",
                                    "drop_sample:",ti$`dropsample (whether to remove the sample for failing quality control)`[r],";",
                                    "race:",ti$race[r],";",
                                    "plate:",ti$`plate (processing plate)`[r],";",
                                    "brain_number:",ti$`brnum (brain number, indicates unique subjects)`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE104210
i=11
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "NA", # Sample type
                             ti$diagnosis[r], # Disease state
                             ti$Sex[r], # Gender
                             ti$age[r], # Age
                             "NA", # Anatomic Location
                             paste0("institute:",ti$institute[r])),# misc
                           nrow=1))
  message(r)
}

# GSE87650
i=12
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`cell type`[r], # Sample type
                             ti$full_diagnosis[r], # Disease state
                             ti$Sex[r], # Gender
                             ti$age[r], # Age
                             "NA", # Anatomic Location
                             paste0("patient_number:",ti$patient_number[r],";",
                                    "age_at_sample:",ti$`age at sample`[r],";",
                                    "smoking_status:",ti$`smoking status`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE111629
i=13
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             ti$`disease state`[r], # Disease state
                             ti$gender[r], # Gender
                             ti$age[r], # Age
                             "NA", # Anatomic Location
                             paste0("ethnicity:",ti$ethnicity[r])),# misc
                           nrow=1))
  message(r)
}

# GSE121633
i=14
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "NA", # Sample type
                             "NA", # Disease state
                             ti$Sex[r], # Gender
                             ti$age[r], # Age
                             "NA", # Anatomic Location
                             "NA"),# misc
                           nrow=1))
  message(r)
}

# GSE93646
i=15
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`sample type`[r], # Sample type
                             "cancer", # Disease state
                             ti$gender[r], # Gender
                             ti$`age (yrs)`[r], # Age
                             "NA", # Anatomic Location
                             paste0("age_info:years;",
                                    "pathology:",ti$pathology[r])),# misc
                           nrow=1))
  message(r)
}

# GSE87648
i=16
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`cell type`[r], # Sample type
                             ti$simplified_diagnosis[r], # Disease state
                             ti$Sex[r], # Gender
                             ti$age[r], # Age
                             "NA", # Anatomic Location
                             paste0("smoking_status:",ti$`smoking status`[r],";",
                                    "patient_number:",ti$patient_number[r],";",
                                    "scan_date:",ti$scan_date[r])),# misc
                           nrow=1))
  message(r)
}

# GSE89278
i=17
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "NA", # Disease state
                             ti$gender[r], # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             "NA"),# misc
                           nrow=1))
  message(r)
}

# GSE71678
i=18
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "NA", # Sample type
                             "NA", # Disease state
                             ti$`infant gender`[r], # Gender
                             ti$`gestational age`[r], # Age
                             "NA", # Anatomic Location
                             paste0("age_info:gestational_age;",
                                    "birth_weight:",ti$`birth weight (grams)`[r],";",
                                    "maternal_age:",ti$`maternal age`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE75248
i=19
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "NA", # Sample type
                             "NA", # Disease state
                             ti$gender[r], # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("batch:",ti$batch[r],";",
                                    "birth_weight_group:",ti$`birth weight group`[r],";",
                                    "attention:",ti$attention[r],";",
                                    "arousal:",ti$arousal[r],";",
                                    "quality_of_movement:",ti$`quality of movement`[r],";",
                                    "lethargy:",ti$lethargy[r])),# misc
                           nrow=1))
  message(r)
}

# GSE84207
i=20
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "NA", # Sample type
                             "NA", # Disease state
                             ti$gender[r], # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("cohort:",ti$cohort[r],";",
                                    "er_status:",ti$`er status`[r],";",
                                    "pam50:",ti$pam50[r])),# misc
                           nrow=1))
  message(r)
}

# GSE51057
i=21
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "NA", # Sample type
                             "NA", # Disease state
                             ti$gender[r], # Gender
                             ti$age[r], # Age
                             "NA", # Anatomic Location
                             paste0("cancer_type:",ti$`cancer type (icd-10)`[r],";",
                                    "age_at_menarche:",ti$`age at menarche`[r],";",
                                    "time_to_diagnosis:",ti$`time to diagnosis`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE61496
i=22
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
sexvar <- ifelse(ti$`sex, 1=m, 2=f`==" 1","male",ifelse(ti$`sex, 1=m, 2=f`==" 2","female","NA"))
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "NA", # Disease state
                             sexvar[r], # Gender
                             ti$age[r], # Age
                             "NA", # Anatomic Location
                             paste0("sex_info:converted_from_numeric;",
                                    "birth_weight:",ti$`birth-weight`[r],";",
                                    "pair_id:",ti$`pair id`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE84493
i=23
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "cancer", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             "NA"),# misc
                           nrow=1))
  message(r)
}

# GSE68838
i=24
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
samplecode <- substr(ti$sample,15,16)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "NA", # Sample type
                             "cancer", # Disease state
                             ti$Sex[r], # Gender
                             ti$age[r], # Age
                             ti$`organism part`[r], # Anatomic Location
                             paste0("sample_type_code:",samplecode[r],";",
                                    "sample_barcode:",ti$sample[r],";",
                                    "colon_polyps_history:",ti$`colon polyps history`[r],";",
                                    "path_m_stage:",ti$`clinical m staging`[r],";",
                                    "path_t_stage:",ti$`pathologic t staging`[r],";",
                                    "path_n_stage:",ti$`pathologic n staging`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE59065
i=25
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`cell/tissue type`[r], # Sample type
                             "NA", # Disease state
                             ti$gender[r], # Gender
                             ti$`age (yr)`[r], # Age
                             "NA", # Anatomic Location
                             paste0("age_info:years;",
                                    "age_group:",ti$agegroup[r])),# misc
                           nrow=1))
  message(r)
}

# GSE72308
i=26
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "cancer", # Disease state
                             "NA", # Gender
                             ti$age_diagnosis[r], # Age
                             "breast", # Anatomic Location
                             paste0("age_info:age_at_diagnosis;",
                                    "cohort_id:",ti$`cohort id`[r],";",
                                    "grade:",ti$grade[r],";",
                                    "subtype_ihc:",ti$subtype_ihc[r],";",
                                    "subtype_pam50:",ti$subtype_pam50[r],";",
                                    "nodal_status:",ti$nodal_status[r],";",
                                    "er_ihc_status:",ti$er_ihc_status[r],";",
                                    "her2_status:",ti$her2_status[r],";",
                                    "relapse_5years:",ti$relapse_5years[r],";",
                                    "date_death:",ti$date_death[r],";",
                                    "cause_death:",ti$cause_death[r])),# misc
                           nrow=1))
  message(r)
}

# GSE60185
i=27
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`sample tissue`[r], # Sample type
                             ti$`disease state`[r], # Disease state
                             ti$gender[r], # Gender
                             ti$age[r], # Age
                             "NA", # Anatomic Location
                             "NA"),# misc
                           nrow=1))
  message(r)
}

# GSE61454
i=28
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             ti$`subject status`[r], # Disease state
                             ti$gender[r], # Gender
                             ti$`subject age`[r], # Age
                             "NA", # Anatomic Location
                             paste0("bmi:",ti$bmi[r])),# misc
                           nrow=1))
  message(r)
}

# GSE101764
i=29
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             ti$`disease state`[r], # Disease state
                             ti$gender[r], # Gender
                             ti$age[r], # Age
                             "NA", # Anatomic Location
                             paste0("patient_id:",ti$patientid[r])),# misc
                           nrow=1))
  message(r)
}

# GSE111223
i=30
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             ti$`disease state`[r], # Disease state
                             ti$gender[r], # Gender
                             ti$age[r], # Age
                             "NA", # Anatomic Location
                             paste0("ethnicity:",ti$ethnicity[r])),# misc
                           nrow=1))
  message(r)
}

# GSE99863
i=31
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "NA", # Disease state
                             ti$Sex[r], # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("season:",ti$season[r])),# misc
                           nrow=1))
  message(r)
}

# GSE85210
i=32
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("smoking_status:",ti$`subject status`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE72872
i=33
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             "NA"),# misc
                           nrow=1))
  message(r)
}

# GSE72874
i=34
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             "NA"),# misc
                           nrow=1))
  message(r)
}

# GSE87640
i=35
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`cell type`[r], # Sample type
                             ti$full_diagnosis[r], # Disease state
                             ti$Sex[r], # Gender
                             ti$`age at sample`[r], # Age
                             "NA", # Anatomic Location
                             paste0("simplified_diagnosis:",ti$simplified_diagnosis[r],";",
                                    "age_info:age_at_sample;",
                                    "patient_number:",ti$patient_number[r])),# misc
                           nrow=1))
  message(r)
}

# GSE97362
i=36
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             ti$`disease state`[r], # Disease state
                             ti$gender[r], # Gender
                             ti$`age (years)`[r], # Age
                             "NA", # Anatomic Location
                             paste0("age_info:years;",
                                    "sample_type:",ti$`sample type`[r],";")),# misc
                           nrow=1))
  message(r)
}

# GSE72021
i=37
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$histology[r], # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("stage:",ti$Stage[r],";",
                                    "grade:",ti$grade[r])),# misc
                           nrow=1))
  message(r)
}

# GSE80261
i=38
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`cell type`[r], # Sample type
                             ti$`disease status`[r], # Disease state
                             ti$Sex[r], # Gender
                             ti$age[r], # Age
                             "NA", # Anatomic Location
                             paste0("self_reported_ethnicity:",ti$`self-reported ethnicity`[r],";",
                                    "disease_subgroup:",ti$`disease subgroup`[r],";",
                                    "date_run:",ti$`date run`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE49149
i=39
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "NA", # Disease state
                             ti$gender[r], # Gender
                             ti$age[r], # Age
                             ti$tissue[r], # Anatomic Location
                             paste0("cellularity:",ti$cellularity[r])),# misc
                           nrow=1))
  message(r)
}

# GSE103186
i=40
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             "NA"),# misc
                           nrow=1))
  message(r)
}

# GSE66351
i=41
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`cell type`[r], # Sample type
                             ti$diagnosis[r], # Disease state
                             ti$Sex[r], # Gender
                             ti$age[r], # Age
                             "brain", # Anatomic Location
                             paste0("brain_region:",ti$brain_region[r],";",
                                    "braak_stage:",ti$braak_stage[r],";",
                                    "braak_region:",ti$brain_region[r],";",
                                    "donor_id:",ti$donor_id[r])),# misc
                           nrow=1))
  message(r)
}

# GSE69138
i=42
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`sample type`[r], # Sample type
                             ti$`disease state`[r], # Disease state
                             ti$gender[r], # Gender
                             ti$age[r], # Age
                             "NA", # Anatomic Location
                             paste0("stroke_subtype:",ti$`stroke subtype`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE66836
i=43
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "NA", # Disease state
                             ti$gender[r], # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("stage:",ti$Stage[r],";",
                                    "p53_status:",ti$`p53 status`[r],";",
                                    "egfr_status:",ti$`egfr status`[r],";",
                                    "kras_status:",ti$`kras status`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE103659
i=44
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "cancer", # Disease state
                             "NA", # Gender
                             ti$age[r], # Age
                             "NA", # Anatomic Location
                             paste0("histology:",ti$histology[r],";",
                                    "mgmt_status:",ti$mgmtstatus[r])),# misc
                           nrow=1))
  message(r)
}

# GSE108462
i=45
ti <- as.data.frame(tls[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "prostate cancer", # Disease state
                             "NA", # Gender
                             ti$age[r], # Age
                             "NA", # Anatomic Location
                             paste0("individual_id:",ti$`individual id`[r],";",
                                    "cohort:",ti$cohort[r],";",
                                    "center:",ti$center[r],";",
                                    "visit:",ti$visit[r],";",
                                    "phenotype:",ti$phenotype[r],";",
                                    "psa_amount:",ti$`psa amount`[r],";",
                                    "treatment_response:",ti$`treatment respsonse`[r])),# misc
                           nrow=1))
  message(r)
}

save(gat,file="geoanno_54gse.rda")
write.csv(gat,file="geoanno_54gse.csv")

#-------
# rep 2
#-------
# filter remaining gse ids and randomize on size
# note: expect larger studies to be clinical populations, 
# smaller studies potentially contain more normal tissue, greater tissue variety

length(unique(gat[,2]))

tlsf <- tls[!names(tls) %in% unique(gat[,2])] # subset
length(tlsf)
tlsf <- tlsf[sample(length(tlsf),length(tlsf))] # randomize order

# GSE62727
i=1
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             ti$condition[r], # Disease state
                             ti$tag[r], # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             "NA"),# misc
                           nrow=1))
  message(r)
}

# GSE70478
i=2
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`cell type`[r], # Sample type
                             "NA", # Disease state
                             ti$gender[r], # Gender
                             ti$age[r], # Age
                             "NA", # Anatomic Location
                             paste0("study_group:",ti$`study group`[r],";",
                                    "matched_pair_id:",ti$`matched pair id`[r],";",
                                    "subject_id:",ti$`subject id`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE70478
i=3
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("cell_line:",ti$`cell line`[r],";",
                                    "mutant:",ti$mutant[r],";",
                                    "treatment:",ti$treatment[r],";",
                                    "duration:",ti$duration[r])),# misc
                           nrow=1))
  message(r)
}

# GSE77955
i=4
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`tissue type`[r], # Sample type
                             "NA", # Disease state
                             ti$gender[r], # Gender
                             ti$`age (yrs)`[r], # Age
                             ti$site[r], # Anatomic Location
                             paste0("age_info:years;",
                                    "tumor_stage:",ti$`tumor stage`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE113775
i=5
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`disease state`[r], # Sample type
                             "NA", # Disease state
                             ti$gender[r], # Gender
                             "NA", # Age
                             ti$tissue[r], # Anatomic Location
                             paste0("study_info:check_patient_disease_state")),# misc
                           nrow=1))
  message(r)
}

# 
i=6
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`sample type`[r], # Sample type
                             ti$`disease state`[r], # Disease state
                             ti$gender[r], # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("cell_line:",ti$`cell line`[r])),# misc
                           nrow=1))
  message(r)
}

# 
i=7
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             "NA"),# misc
                           nrow=1))
  message(r)
}

# 
i=8
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "NA", # Sample type
                             "NA", # Disease state
                             ti$Sex[r], # Gender
                             ti$age[r], # Age
                             "NA", # Anatomic Location
                             paste0("race:",ti$race[r],";",
                                    "ischemia:",ti$ischemia[r])),# misc
                           nrow=1))
  message(r)
}

# 
i=9
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`cell type`[r], # Sample type
                             ti$`disease status`[r], # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("expansion_step:",ti$`expansion step`[r])),# misc
                           nrow=1))
  message(r)
}

# 
i=10
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`cell type`[r], # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("study_info:enriched_primary_cells",
                                    "cell_line:",ti$`cell line`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE67349
i=11
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "NA", # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("study_info:cell_line_study",
                                    "cell_line:",ti$`cell line`[r],";",
                                    "treatment:",ti$treatment[r])),# misc
                           nrow=1))
  message(r)
}

# study 66: 
i=12
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "cancer", # Disease state
                             ti$gender[r], # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("patient_id:",ti$`patient id`[r],";",
                                    "year_tissue_collection:",ti$`year of tissue collection`[r],";",
                                    "preservation_method:",ti$`specimen preservation method`[r],";",
                                    "input_amount:",ti$`input amount`[r])),# misc
                           nrow=1))
  message(r)
}

save(gat,file="geoanno_66gse.rda")
write.csv(gat,file="geoanno_66gse.csv")

#------
# rep3
#------
length(unique(gat[,2])) # 66
tlsf <- tls[!names(tls) %in% unique(gat[,2])] # subset
length(tlsf) # 412
tlsf <- tlsf[sample(length(tlsf),length(tlsf))] # randomize order

# GSE92909
i=1
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "cell_line", # Sample type
                             ti$`disease subtype`[r], # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("study_info:cell_line_study;",
                                    "cell_line:",ti$`cell line`[r],";",
                                    "treatment:",ti$treatment[r])),# misc
                           nrow=1))
  message(r)
}

# GSE93266 
i=2
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "NA", # Sample type
                             ti$`disease state`[r], # Disease state
                             ti$Sex[r], # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             "NA"),# misc
                           nrow=1))
  message(r)
}

# GSE75434 
i=3
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("individual:",ti$individual[r])),# misc
                           nrow=1))
  message(r)
}

# GSE103768 
i=4
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "NA", # Disease state
                             ti$gender[r], # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("time:",ti$time[r],";",
                                    "outcome:",ti$outcome[r])),# misc
                           nrow=1))
  message(r)
}

# GSE66564  
i=5
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`cell type`[r], # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             ti$tissue[r], # Anatomic Location
                             paste0("donor_id:",ti$`donor id`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE65638  
i=6
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             ti$age[r], # Age
                             ti$tissue[r], # Anatomic Location
                             paste0("pair_id_number:",ti$`pair id number`[r],";",
                                    "ethnicity:",ti$ethnicity[r])),# misc
                           nrow=1))
  message(r)
}

# GSE53261  
i=7
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`cell type`[r], # Sample type
                             "NA", # Disease state
                             ti$gender[r], # Gender
                             "NA", # Age
                             ti$tissue[r], # Anatomic Location
                             paste0("cell_line:",ti$`cell line`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE67485   
i=8
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "NA", # Disease state
                             ti$gender[r], # Gender
                             ti$age[r], # Age
                             ti$tissue[r], # Anatomic Location
                             paste0("cell_line:",ti$`cell line`[r],";",
                                    "agent:",ti$agent[r])),# misc
                           nrow=1))
  message(r)
}

# GSE87424  
i=9
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`cell line`[r], # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("treatment:",ti$treatment[r],";",
                                    "control_vendor:",ti$`control vendor`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE78277   
i=10
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`cell type`[r], # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "foreskin", # Anatomic Location
                             paste0("transformation_stage:",ti$`transformation stage`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE114935   
i=11
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "NA", # Sample type
                             "NA", # Disease state
                             ti$characteristics[r], # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("tp:",ti$tp[r],";",
                                    "cd8t_adult:",ti$cd8t_adult[r],";",
                                    "cd4t_adult:",ti$cd4t_adult[r],";",
                                    "nk_adult:",ti$nk_adult[r],";",
                                    "bcell_adult:",ti$bcell_adult[r],";",
                                    "maternal_age:",ti$maternalage[r])),# misc
                           nrow=1))
  message(r)
}

# GSE81228   
i=12
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
ttypevar <- ifelse(ti$`tumor type`=="NA",ti$`cell type`,ti$`tumor type`)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ttypevar[r], # Sample type
                             "NA", # Disease state
                             ti$gender[r], # Gender
                             ti$`patient age`[r], # Age
                             "NA", # Anatomic Location
                             paste0("tumor_histology:",ti$`tumor histology`[r],";",
                                    "stage:",ti$Stage[r],";",
                                    "grade:",ti$grade[r])),# misc
                           nrow=1))
  message(r)
}

# GSE61441   
i=13
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`sample type`[r], # Sample type
                             ti$`disease state`[r], # Disease state
                             ti$gender[r], # Gender
                             "NA", # Age
                             ti$tissue[r], # Anatomic Location
                             paste0("NA")),# misc
                           nrow=1))
  message(r)
}

# GSE73901   
i=14
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`cell type`[r], # Sample type
                             "NA", # Disease state
                             ti$gender[r], # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("population:",ti$population[r])),# misc
                           nrow=1))
  message(r)
}

# GSE105260  
i=15
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("NA")),# misc
                           nrow=1))
  message(r)
}

# GSE80226  
i=16
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`cell type`[r], # Sample type
                             "NA", # Disease state
                             ti$gender[r], # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("cell_line:",ti$`cell line`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE72867   
i=17
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`cell type`[r], # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             ti$`tissue type`[r], # Anatomic Location
                             paste0("NA")),# misc
                           nrow=1))
  message(r)
}

# GSE54375    
i=18
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`cell type`[r], # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "lung", # Anatomic Location
                             paste0("cell_line:",ti$`cell line`[r],";",
                                    "treatment:",ti$`treated with`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE109430    
i=19
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`cell type`[r], # Sample type
                             ti$`disease state`[r], # Disease state
                             ti$gender[r], # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("NA")),# misc
                           nrow=1))
  message(r)
}

# GSE117439   
i=20
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             ti$`disease state`[r], # Disease state
                             ti$Sex[r], # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("group:",ti$group[r],";",
                                    "er_status:",ti$`estrogen receptor status`[r],";",
                                    "tumor_status:",ti$`tumor status`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE80369   
i=21
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "NA", # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("culture_methods:",ti$`culture methods`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE66077   
i=22
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`cell line`[r], # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("stage:",ti$Stage[r],";",
                                    "group:",ti$group[r])),# misc
                           nrow=1))
  message(r)
}

# GSE107460   
i=23
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "NA", # Sample type
                             ti$`sample type`[r], # Disease state
                             ti$Sex[r], # Gender
                             ti$age[r], # Age
                             "blood", # Anatomic Location
                             paste0("race:",ti$race[r],";",
                                    "fetal_intolerance:",ti$`fetal intolerance`[r],";",
                                    "visit:",ti$visit[r],";",
                                    "ga_sample:",ti$`ga sample`[r],";",
                                    "pregnancy_complication:",ti$`pregnancy complication`[r])),# misc
                           nrow=1))
  message(r)
}

# GSE95816   
i=24
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("predicted_gender:",ti$predicted_gender[r],";",
                                    "dna_conc:",ti$dna_concentration[r],";",
                                    "volume:",ti$volume[r],";",
                                    "input:",ti$input[r])),# misc
                           nrow=1))
  message(r)
}

# GSE77201  
i=25
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "cell_line", # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("cell_line:",ti$`cell line`[r],";",
                                    "transfection:",ti$transfection[r])),# misc
                           nrow=1))
  message(r)
}

# GSE92577    
i=26
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             "cancer", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "brain", # Anatomic Location
                             paste0("NA")),# misc
                           nrow=1))
  message(r)
}

# GSE107351   
i=27
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "NA", # Sample type
                             "NA", # Disease state
                             ti$gender[r], # Gender
                             "NA", # Age
                             "blood", # Anatomic Location
                             paste0("tissue_info:",ti$tissue[r])),# misc
                           nrow=1))
  message(r)
}

# GSE55491  
i=28
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             ti$`disease state`[r], # Disease state
                             ti$gender[r], # Gender
                             "NA", # Age
                             "blood", # Anatomic Location
                             paste0("NA")),# misc
                           nrow=1))
  message(r)
}

# GSE81794   
i=29
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             "NA", # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("treatment:",ti$treatment[r])),# misc
                           nrow=1))
  message(r)
}

# GSE56596  
i=30
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`cell type`[r], # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("tissue_type:",ti$`tissue type`[r],";",
                                    "subtype:",ti$subtype[r])),# misc
                           nrow=1))
  message(r)
}

# GSE74432   
i=31
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             ti$`disease state`[r], # Disease state
                             ti$gender[r], # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("NA")),# misc
                           nrow=1))
  message(r)
}

# GSE94987   
i=32
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`sample type`[r], # Sample type
                             "NA", # Disease state
                             "NA", # Gender
                             "NA", # Age
                             "NA", # Anatomic Location
                             paste0("NA")),# misc
                           nrow=1))
  message(r)
}

# GSE71955   
i=33
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$`cell type`[r], # Sample type
                             ti$diagnosis[r], # Disease state
                             ti$gender[r], # Gender
                             ti$age[r], # Age
                             ti$tissue[r], # Anatomic Location
                             paste0("NA")),# misc
                           nrow=1))
  message(r)
}

# GSE77276  
i=34
ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
head(ti)
for(r in 1:nrow(ti)){
  gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                             ti$gse[r], # GSE ID
                             ti$tissue[r], # Sample type
                             ti$`subject status`[r], # Disease state
                             ti$gender[r], # Gender
                             ti$age[r], # Age
                             "NA", # Anatomic Location
                             paste0("patient_id:",ti$`patient id`[r],";",
                                    "tumor_size:",ti$`tumor size (cm)`[r],";",
                                    "hbv_infection:",ti$`hbv infection`[r])),# misc
                           nrow=1))
  message(r)
}

save(gat,file="geoanno_100gse.rda")
write.csv(gat,file="geoanno_100gse.csv")

#-------
# rep 4
#-------
length(unique(gat[,2])) # 100
tlsf <- tls[!names(tls) %in% unique(gat[,2])] # subset
length(tlsf) # 378
tlsf <- tlsf[sample(length(tlsf),length(tlsf))] # randomize order

# study 101:GSE100386   
{
  i=1
  ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("atopic:",ti$atopic[r])),# misc
                             nrow=1))
    message(r)
  }
}

# study 102:GSE102177   
{
  i=2
  ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("group:",ti$group[r],";",
                                      "sibling:",ti$sibling[r])),# misc
                             nrow=1))
    message(r)
  }
}

# study 103:GSE102996   
{
  i=3
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# study 104: GSE96866   
{
  i=4
  ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell type`[r],";",
                                      "transfected_with:",ti$`transfected with`[r],";",
                                      "time_point:",ti$`time point`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# study 105:GSE87053   
{
  i=5
  ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`tissue type`[r], # Sample type
                               "cancer", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("hpv_status:",ti$`hpv status`[r],";",
                                      "subject_status:",ti$`subject status`[r],";",
                                      "study_group:",ti$`disease state`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# study 106: GSE65307   
{
  i=6
  ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("diagnosis:",ti$diagnosis[r])),# misc
                             nrow=1))
    message(r)
  }
}

# study 107: GSE86648  
{
  i=7
  ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# study 108: GSE94368   
{
  i=8
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`tissue type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# study 109: GSE94462  
{
  i=9
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$`disease status`[r], # Disease state
                               ti$Sex[r], # Gender
                               ti$age[r], # Age
                               "cornea", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# study 110: GSE81438   
{
  i=10
  ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  ttvar <- ifelse(ti$tumor=="NA",ti$`cell line`,ti$`tumor type`)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ttvar[r], # Sample type
                               "cancer", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("tumor_type:",ti$`tumor type`[r],";",
                                      "cell_line:",ti$`cell line`[r],";",
                                      "tumor_generation:",ti$`tumor generation`[r],";",
                                      "treatment:",ti$treatment[r],";",
                                      "pdx_model:",ti$`pdx model`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# study 111:GSE111396    
{
  i=11
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$disease[r], # Disease state
                               ti$Sex[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("passage:",ti$passage[r],";",
                                      "smoker:",ti$smoker[r],";",
                                      "pack_years:",ti$`pack years`[r],";",
                                      "gold_stage",ti$`gold stage`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# study 112: GSE71626  
{
  i=12
  ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample type`[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# study 113: GSE71458
{
  i=13
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "differentiation_status:",ti$`differentiation status`[r],";",
                                      "passage:",ti$passage[r],";",
                                      "cell_cycle_phase:",ti$`cell cycle phase`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# study 114:GSE81215
{
  i=14
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# study 115: GSE94326
{
  i=15
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`tissue type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# study 116:GSE65467 
{
  i=16
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("treatment:",ti$treatment[r],";",
                                      "passage:",ti$passage[r])),# misc
                             nrow=1))
    message(r)
  }
}

# study 117:GSE83379
{
  i=17
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "time:",ti$time[r],";",
                                      "infection:",ti$infection[r])),# misc
                             nrow=1))
    message(r)
  }
}

# study 118:GSE85938
{
  i=18
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("passage:",ti$passage[r],";",
                                      "dox:",ti$dox[r],";",
                                      "genotype:",ti$genotype[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 119:GSE79185
{
  i=19
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               ti$tissue[r], # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 120:GSE109446
{
  i=20
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$diagnosis[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "nasal epithelium", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 121:GSE88940
{
  i=21
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("ethnicity:",ti$ethnicity[r],";",
                                      "bmi:",ti$bmi[r],";",
                                      "weight_kg:",ti$`weight (kg)`[r],";",
                                      "height_cm:",ti$`height (cm)`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 122:GSE80243
{
  i=22
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "agent:",ti$agent[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 123:GSE85087
{
  i=23
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  ttvar <- ifelse(ti$`tissue type (primary=0, recurrent=1, recurrent2=2)`==" 0","primary",
                  ifelse(ti$`tissue type (primary=0, recurrent=1, recurrent2=2)`==" 1","recurrent1",
                         ifelse(ti$`tissue type (primary=0, recurrent=1, recurrent2=2)`==" 2","recurrent2","NA")))
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ttvar[r], # Sample type
                               "cancer", # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("survival:",ti$`survival (months`[r],";",
                                      "time_to_relapse_mo:",ti$`time to relapse (months)`[r],";",
                                      "treatment_prior_to_second_surgery:",ti$`treatment prior to second surgery`[r],";",
                                      "replicate_array:",ti$`replicate array (0=no, 1=yes)`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 124:GSE61279
{
  i=124
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  dxvar <- ifelse(ti$group==" Cancer","cancer","NA")
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "liver", # Sample type
                               dxvar[r], # Disease state
                               ti$Sex[r], # Gender
                               ti$agegroup[r], # Age
                               "liver", # Anatomic Location
                               paste0("age_info:age_group;",
                                      "group:",ti$group[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 125:GSE99553
{
  i=125
  ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  dxvar <- ifelse(ti$`cancer status`==" Case","cancer","control")
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               dxvar[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("h_pylori_infection:",ti$`h. pylori infection`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 126:GSE81211
{
  i=26
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";")),# misc
                             nrow=1))
    message(r)
  }
}

# 127:GSE77348
{
  i=127
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 128:GSE89474
{
  i=28
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               ti$`disease status`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$`age (years)`[r], # Age
                               "blood", # Anatomic Location
                               paste0("age_info:years;",
                                      "twin_pair:",ti$`twin pair`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 129: GSE80468
{
  i=29
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$diagnosis[r], # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "blood", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 130: GSE61446
{
  i=130
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "liver", # Sample type
                               ti$`subject status`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$`subject age`[r], # Age
                               "liver", # Anatomic Location
                               paste0("bmi:",ti$bmi[r],";")),# misc
                             nrow=1))
    message(r)
  }
}

# 131: GSE80559
{
  i=131
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("sleep:",ti$sleep[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 132: GSE80969
{
  i=132
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 133:GSE85467
{
  i=133
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`tissue subtype`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "gastric_tissue", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 134:GSE100134
{
  i=134
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 135: GSE100561
{
  i=135
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell population`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "blood", # Anatomic Location
                               paste0("ethnic_group:",ti$`ethnic group`[r],";",
                                      "individual:",ti$individual[r],";",
                                      "infection:",ti$infection[r],";",
                                      "cell_type:",ti$`cell type`[r]
                                      )),# misc
                             nrow=1))
    message(r)
  }
}

# 136:GSE72245
{
  i=136
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "breast", # Anatomic Location
                               paste0("age_info:age_diagnosis;",
                                      "cohort_id:",ti$`cohort id`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 137:GSE71523
{
  i=137
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample type`[r], # Sample type
                               "NA", # Disease state
                               ti$tag[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("source_id:",ti$`sample source id`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 138:GSE83933
{
  i=138
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "cancer", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("disease_state:",ti$`disease state`[r],";",
                                      "who_grade:",ti$`who grade`[r],";",
                                      "history:",ti$history[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 139:GSE112047
{
  i=139
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "prostate cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("source_id:",ti$`source id`[r],";",
                                      "tumor_cellularity:",ti$`tumor cellularity`[r],";",
                                      "sample_preparation:",ti$`sample preparation`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 140:GSE71457
{
  i=140
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "differentiation_status:",ti$`differentiation status`[r],";",
                                      "cell_cycle_phase:",ti$`cell cycle phase`[r],";",
                                      "passage:",ti$passage[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 141:GSE68344
{
  i=141
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 142:GSE80368
{
  i=142
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("culture_methods:",ti$`culture methods`[r])),# misc
                             nrow=1))
    message(r)
  }
}


save(gat,file="geoanno_142gse.rda")
write.csv(gat,file="geoanno_142gse.csv")

#-------
# rep 5
#-------

length(unique(gat[,2])) # 142
tlsf <- tls[!names(tls) %in% unique(gat[,2])] # subset
length(tlsf) # 336
tlsf <- tlsf[sample(length(tlsf),length(tlsf))] # randomize order

# 143:GSE72354
{
  i=1
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 144:GSE102119
{
  i=2
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "ovary", # Anatomic Location
                               paste0("patient_id:",ti$`patient id`[r],";",
                                      "treatment:",ti$treatment[r],";")),# misc
                             nrow=1))
    message(r)
  }
}

# 145:GSE100197
{
  i=3
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample tissue`[r], # Sample type
                               "NA", # Disease state
                               ti$`fetal sex`[r], # Gender
                               ti$`gestational age`[r], # Age
                               "NA", # Anatomic Location
                               paste0("subject_id:",ti$`subject id`[r],";",
                                      "pathology_group:",ti$`pathology group`[r],";",
                                      "450k_plate:",ti$`450k plate`[r],";",
                                      "age_info:gestational_age")),# misc
                             nrow=1))
    message(r)
  }
}

# 146:GSE61453
{
  i=146
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$`subject status`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$`subject age`[r], # Age
                               "NA", # Anatomic Location
                               paste0("bmi:",ti$bmi[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 147:GSE61160
{
  i=4
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("diagnosis:",ti$diagnosis[r],";",
                                      "idh_status:",ti$`idh status`[r])),# misc
                             nrow=1))
    message(r)
  }
}

#148: GSE89778
{
  i=5
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample  group`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("subject:",ti$subject[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 149: GSE85529
{
  i=6
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 150: GSE81309
{
  i=7
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 151:GSE65163
{
  i=8
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  racevar <- ifelse(ti$race_a_carrib==" Yes","race_a_carrib",
                    ifelse(ti$race_aa==" Yes","race_aa",
                           ifelse(ti$race_asian==" Yes","race_asian",
                                  ifelse(ti$race_hispanic==" Yes","race_hispanic",
                                         ifelse(ti$race_indian_alaska==" Yes","race_indian_alaska",
                                                ifelse(ti$race_white==" Yes","race_white","NA"))))))
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "nasal epithelium", # Anatomic Location
                               paste0("site:",ti$site[r],";",
                                      "race:",racevar[r],";",
                                      "asthma:",ti$asthma[r],";",
                                      "lnige:",ti$lnige[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 152: GSE71328
{
  i=152
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               ti$tissue[r], # Anatomic Location
                               paste0("individual:",ti$individual[r],";",
                                      "dna_conversion_protocol:",ti$`dna conversion protocol`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 153: GSE73549
{
  i=9
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`tissue type`[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               ti$`tissue type`[r], # Anatomic Location
                               paste0("patient_id:",ti$`patient id`[r],";",
                                      "sample_type:",ti$`sample type`[r],";",
                                      "sample_info:",ti$`sample info`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 154: GSE84397
{
  i=10
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "colon", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "transfected_with:",ti$`transfected with`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 155: GSE61452
{
  i=11
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$`subject status`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$`subject age`[r], # Age
                               "NA", # Anatomic Location
                               paste0("bmi:",ti$bmi[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 156:GSE89251
{
  i=12
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               ti$`age (ys)`[r], # Age
                               "NA", # Anatomic Location
                               paste0("age_info:years;",
                                      "clinical_activity:",ti$`clinical activity`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 157: GSE67733
{
  i=13
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 158: GSE102994
{
  i=14
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 159: GSE57992
{
  i=15
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$pluripotent[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 160: GSE122126
{
  i=16
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample type`[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 161: GSE109096
{
  i=17
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("ischemia:",ti$ischemia[r],";",
                                      "race:",ti$race[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 162: GSE109914
{
  i=18
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               ti$`disease status`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("arsenic_exposure:",ti$`arsenic exposure`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 163: GSE71825
{
  i=163
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "blood", # Anatomic Location
                               paste0("activation_state:",ti$`activation state`[r],";",
                                      "activation_protocol:",ti$`activation protocol`[r],";",
                                      "cell_source:",ti$`cell source`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 164:GSE116924
{
  i=19
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("vector_type:",ti$`vector type`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 165: GSE66078
{
  i=20
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "group:",ti$group[r],";",
                                      "stage:",ti$Stage[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 166: GSE73515
{
  i=21
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               ti$`age at diagnosis`[r], # Age
                               "NA", # Anatomic Location
                               paste0("age_info:age_at_diagnosis",
                                      "inss_stage:",ti$`inss stage`[r],";",
                                      "current_risk_category:",ti$`current risk category`[r],";",
                                      "mycn_status:",ti$`mycn status`[r],";",
                                      "1p_status:",ti$`1p status`[r],";",
                                      "11q_status:",ti$`11q status`[r],";",
                                      "17q_status:",ti$`17q status`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 167:GSE84400
{
  i=167
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "colon", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "transfected_with:",ti$`transfected with`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 168: GSE62929
{
  i=22
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               ti$`cell line source tissue`[r], # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                 "genotype_variation:",ti$`genotype/variation`[r],";",
                                      "cell_line_background:",ti$`cell line background`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 169: GSE85688
{
  i=23
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "colon", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 170: GSE97529
{
  i=24
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$diagnosis[r], # Disease state
                               ti$Sex[r], # Gender
                               "NA", # Age
                               "bone", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 171: GSE112696
{
  i=25
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$diagnosis[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               ti$tissue[r], # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 172: GSE108567
{
  i=26
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$`fetal sex`[r], # Gender
                               ti$`gestational age (week, rounded)`[r], # Age
                               "placenta", # Anatomic Location
                               paste0("age_info:gestational_age_week_rounded;",
                                      "gender_info:fetal_gender;",
                                      "processing_group:",ti$`processing group`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 173: GSE83691
{
  i=27
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 174: GSE72364
{
  i=28
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 175: GSE66295
{
  i=29
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 176: GSE77954
{
  i=30
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`tissue type`[r], # Sample type
                               "cancer", # Disease state
                               ti$gender[r], # Gender
                               ti$`age (yrs)`[r], # Age
                               ti$site[r], # Anatomic Location
                               paste0("age_info:years;",
                                      "tumor_stage:",ti$`tumor stage`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 177:GSE103010
{
  i=31
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "bone marrow, mouse", # Anatomic Location
                               paste0("host_mouse_strain:",ti$`host mouse strain`[r],";",
                                      "humanized_duration:",ti$`humanized duration`[r],";",
                                      "primary_cell_type:",ti$`primary cell type`[r],";",
                                      "primary_tissue_type:",ti$`primary tissue type`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 178: GSE81008
{
  i=32
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "genotype:",ti$genotype[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 179: GSE63106
{
  i=33
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               ti$age[r], # Age
                               paste0(ti$joint[r]," cartilage"), # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 180: GSE70737
{
  i=34
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 181: GSE104087
{
  i=35
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$diagnosis[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "nasal epithelium", # Anatomic Location
                               paste0("race:",ti$race[r],";",
                                      "time:",ti$time[r],";",
                                      "length_of_stay:",ti$`length of stay`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 182: GSE83320
{
  i=36
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "skin", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 183: GSE99624
{
  i=37
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "whole blood", # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "whole blood", # Anatomic Location
                               paste0("donor_id:",ti$`donor id`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 184: GSE112877
{
  i=38
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("group:",ti$group[r],";",
                                      "stability:",ti$stability[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 185: GSE107205
{
  i=39
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("diet:",ti$diet[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 186: GSE106089
{
  i=40
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "placenta", # Anatomic Location
                               paste0("maternal_race:",ti$`maternal race`[r],";",
                                      "fetal_sex:",ti$`fetal sex`[r],";",
                                      "antibiotic_use:",ti$`antibiotic use`[r],";",
                                      "vaginal_infection:",ti$`vaginal infection`[r],";",
                                      "uti:",ti$uti[r],";",
                                      "c_section:",ti$`c-section`[r],";",
                                      "bacteria_present:",ti$`bacteria present`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 187: GSE72120
{
  i=41
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "saliva", # Anatomic Location
                               paste0("pregnancy:",ti$pregancy[r],";",
                                      "pma_birth:",ti$pmabirth[r],";",
                                      "steroids:",ti$steroids[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 188: GSE68342
{
  i=42
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 189: GSE92378
{
  i=43
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";", 
                                      "ebv_status:",ti$`ebv status`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 190: GSE98876
{
  i=44
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("sample_group:",ti$`sample group`[r],";",
                                      "cigarettes_per_day:",ti$`cigarettes per day`[r],";",
                                      "nationality:",ti$nationality[r],";",
                                      "days_since_last_drink:",ti$`days since last drink`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 191: GSE77797
{
  i=45
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("study_info:wbc_compositions")),# misc
                             nrow=1))
    message(r)
  }
}

# 192: GSE102468
{
  i=46
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "blood", # Sample type
                               ti$diagnosis[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "blood", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 193: GSE86402
{
  i=47
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell line`[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "colon", # Anatomic Location
                               paste0("genotype_variant:",ti$`genotype/variant`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 194: GSE77135
{
  i=48
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell line`[r], # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("race:",ti$race[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 195: GSE85566
{
  i=49
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$`disease status`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("ethnicity:",ti$ethnicity[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 196: GSE67484
{
  i=50
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               ti$tissue[r], # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 197: GSE75443
{
  i=51
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$disease[r], # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("radiation_response:",ti$`radiation response`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 198: GSE112314
{
  i=52
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "saliva", # Anatomic Location
                               paste0("lsi_score:",ti$`lsi score`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 199: GSE93589
{
  i=53
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`tissue type`[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("subject:",ti$subject[r],";",
                                      "race:",ti$race[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 200: GSE105066
{
  i=54
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "msssi:",ti$msssi[r],";",
                                      "day:",ti$day[r],";",
                                      "treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

save(gat,file="geoanno_200gse.rda")
write.csv(gat,file="geoanno_200gse.csv")

#-------
# rep 6
#-------

length(unique(gat[,2])) # 200
tlsf <- tls[!names(tls) %in% unique(gat[,2])] # subset
length(tlsf) # 278
tlsf <- tlsf[sample(length(tlsf),length(tlsf))] # randomize order

# 201: GSE73115
{
  i=1
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("individual_id:",ti$`individual id`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 202: GSE69502
{
  i=2
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample tissue`[r], # Sample type
                               ti$`ntd status`[r], # Disease state
                               ti$`fetal sex`[r], # Gender
                               ti$`fetal gestational age (weeks)`[r], # Age
                               ti$`sample tissue`[r], # Anatomic Location
                               paste0("ntd_status:",ti$`ntd status`[r],";",
                                      "age_info:fetal_gestational_age_weeks;",
                                      "gender_info:fetal_gender")),# misc
                             nrow=1))
    message(r)
  }
}

# 203:GSE74548
{
  i=3
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$`age at baseline`[r], # Age
                               "NA", # Anatomic Location
                               paste0("intervention:",ti$intervention[r],";",
                                      "age_info:age_at_baseline;",
                                      "subject:",ti$subject[r],";",
                                      "time_point:",ti$`time point`[r],";",
                                      "bisulphite_plate:",ti$`bisulfite plate`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 204: GSE80685
{
  i=4
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "prostate", # Anatomic Location
                               paste0("disease_state:",ti$`disease state`[r],";",
                                      "sample:",ti$`^sample`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 205: GSE107353
{
  i=5
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 206: GSE83917
{
  i=6
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "prostate", # Anatomic Location
                               paste0("disease_state:",ti$`disease state`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 207: GSE84043
{
  i=7
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "prostate", # Anatomic Location
                               paste0("disease_state:",ti$`disease state`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 208:GSE114753 
{
  i=8
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "sperm", # Anatomic Location
                               paste0("sperm_concentration_m_per_ml:",ti$`sperm concentration (m/ml)`[r],";",
                                      "pack_years:",ti$`pack years`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 209: GSE70460
{
  i=9
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$`age at diagnosis`[r], # Age
                               "NA", # Anatomic Location
                               paste0("brain_localization:",ti$`localization within the brain`[r],";",
                                      "age_info:age_at_diagnosis")),# misc
                             nrow=1))
    message(r)
  }
}

# 210:GSE86078 
{
  i=10
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "colorectal", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "supplier:",ti$supplier[r],";",
                                      "catalog_id:",ti$catalog_id[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 211:GSE102120
{
  i=11
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "ovary", # Anatomic Location
                               paste0("patient_id:",ti$`patient id`[r],";",
                                      "treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 212: GSE68825
{
  i=12
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`histologic type`[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$Sex[r], # Gender
                               ti$age[r], # Age
                               ti$`organism part`[r], # Anatomic Location
                               paste0("pack_years_smoked:",ti$`number pack years smoked`[r],";",
                                      "disease_staging:",ti$`disease staging`[r],";",
                                      "smoking_history:",ti$`smoking history`[r],";",
                                      "sample_barcode:",ti$sample[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 213: GSE65183
{
  i=13
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "skin", # Anatomic Location
                               paste0("mapki_sensitivity:",ti$`mapki sensitivity`[r],";",
                                      "mapki_treatment:",ti$`mapki treatment`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 214: GSE65186
{
  i=14
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "skin", # Anatomic Location
                               paste0("mapki_sensitivity:",ti$`mapki sensitivity`[r],";",
                                      "mapki_treatment:",ti$`mapki treatment`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 215: GSE97466
{
  i=15
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$histology[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("race:",ti$race[r],";",
                                      "variant:",ti$variant[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 216: GSE79257
{
  i=16
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "neonatal_blood", # Anatomic Location
                               paste0("conception_type:",ti$`conception type`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 217:GSE89253 
{
  i=17
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               ti$`age (ys)`[r], # Age
                               "NA", # Anatomic Location
                               paste0("age_info:years;",
                                      "clinical_activity:",ti$`clinical activity`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 218: GSE71957
{
  i=18
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$diagnosis[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               ti$tissue[r], # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 219: GSE104293
{
  i=19
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               ti$disease.state[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("patient:",ti$patient[r],";",
                                      "batch:",ti$batch[r],";",
                                      "mgmt_score:", ti$mgmt_score[r],";",
                                      "mgmt_status:",ti$mgmt_status[r],";",
                                      "molecular_subtype:",ti$molecular_subtype[r],";",
                                      "treatment:",ti$treatment[r],";",
                                      "tissue_type:",ti$tissue_type[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 220: GSE107459
{
  i=20
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("race:",ti$race[r],";",
                                      "fetal_intolerance:",ti$`fetal intolerance`[r],";",
                                      "visit:",ti$visit[r],";",
                                      "pregnancy_complication:",ti$`pregnancy complication`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 221: GSE93933
{
  i=21
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$`sample group`[r], # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "blood", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 222:GSE79009
{
  i=22
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               ti$`disease state`[r], # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 223:GSE87095
{
  i=23
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               ti$`disease state`[r], # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 224: GSE65820
{
  i=24
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               ti$tissue[r], # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 225: GSE65821
{
  i=25
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               ti$tissue[r], # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 226: GSE101961
{
  i=26
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("race:",ti$race[r],";",
                                      "bmi:",ti$bmi[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 227: GSE68060
{
  i=27
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("type:",ti$type[r],";",
                                      "country:",ti$country[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 228: GSE72251
{
  i=28
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               ti$age_diagnosis[r], # Age
                               "breast", # Anatomic Location
                               paste0("subtype_ihc:",ti$subtype_ihc[r],";",
                                      "cohort_id:",ti$`cohort id`[r],";",
                                      "age_info:age_diagnosis;",
                                      "grade:",ti$grade[r],";",
                                      "cause_death:",ti$cause_death[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 229: GSE67393
{
  i=29
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "blood", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 230:GSE85568
{
  i=30
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$`disease status`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "airway_epithelial_cells", # Anatomic Location
                               paste0("ethnicity:",ti$ethnicity[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 231:GSE85828
{
  i=31
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("differentiation_stage:",ti$`differentiation stage`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 232: GSE61278
{
  i=32
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               ti$agegroup[r], # Age
                               "liver", # Anatomic Location
                               paste0("group:",ti$group[r],";",
                                      "age_info:age_group")),# misc
                             nrow=1))
    message(r)
  }
}

# 233: GSE105123
{
  i=33
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("subject:",ti$subject[r],";",
                                 "time_point:",ti$`time point`[r],";",
                                      "height_in:",ti$`height (in.)`[r],";",
                                      "weight_kg:",ti$`weight (kg.)`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 234: GSE105124
{
  i=34
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("subject:",ti$subject[r],";",
                                      "time_point:",ti$`time point`[r],";",
                                      "height_in:",ti$`height (in.)`[r],";",
                                      "weight_kg:",ti$`weight (kg.)`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 235: GSE73518
{
  i=35
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               ti$`age at diagnosis`[r], # Age
                               "NA", # Anatomic Location
                               paste0("age_info:age_at_diagnosis;",
                                      "current_risk_category:",ti$`current risk category`[r],";",
                                      "inss_stage:",ti$`inss stage`[r],";",
                                      "mycn_status:",ti$`mycn status`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 236: GSE47915
{
  i=236
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "prostate", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 237: GSE69636
{
  i=37
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("gestational_age:",ti$`gestational age`[r],";",
                                      "socioeconomic_score:",ti$`socioeconomic score`[r],";",
                                      "smoke_ever:",ti$`smoke ever`[r],";",
                                      "birth_weight:",ti$`birth weight`[r],";",
                                      "cell_line:",ti$`cell line`[r],";",
                                      "agent:",ti$agent[r],";",
                                      "dose:",ti$dose[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 238: GSE88883
{
  i=38
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               ti$`subject age`[r], # Age
                               "breast", # Anatomic Location
                               paste0("bmi:",ti$`body mass index`[r],";",
                                      "pregnancy_ever:",ti$`pregnancy (ever)`[r],";",
                                      "family_history_breast_or_ovarian_cancer:",ti$`family history (first-degree relative with breast or ovarian cancer)`[r],";",
                                      "gail_model_risk_score:",ti$`gail risk model score (for women over 36 years old)`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 239: GSE63695
{
  i=39
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("scan_date:",ti$`scan date`[r],";",
                                      "kl_score:",ti$`kl score`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 240: GSE72556
{
  i=40
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("adult_age:",ti$`adult age`[r],";",
                                      "child_age:",ti$`child age`[r],";",
                                      "adult_bmi:",ti$`adult bmi`[r],";",
                                      "child_bmi:",ti$`child bmi`[r],";",
                                      "child_gender:",ti$`child gender`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 241: GSE69633
{
  i=41
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$`gestational age`[r], # Age
                               "umbilical_cord_blood", # Anatomic Location
                               paste0("age_info:gestational_age;",
                                      "smoke_ever:",ti$`smoke ever`[r],";",
                                      "birth_weight:",ti$`birth weight`[r],";",
                                      "socioeconomic_score:",ti$`socioeconomic score`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 242:GSE75153
{
  i=42
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("disease_state:",ti$`disease state`[r],";",
                                      "molecular_subtype:",ti$`molecular subtype`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 243: GSE67170
{
  i=43
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 244:GSE108576
{
  i=44
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("sample_type:",ti$`sample type`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 245:GSE98203
{
  i=45
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("cohort:",ti$cohort[r],";",
                                      "cause_of_death:",ti$`Cause of death`[r],";",
                                      "race:",ti$race[r],";",
                                      "tissue_ph:",ti$tissue_ph[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 246:GSE120250
{
  i=46
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("outlier:",ti$outlier[r],";",
                                      "art_treatment:",ti$`art treatment`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 247:GSE79556
{
  i=47
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "tongue_carcinoma", # Anatomic Location
                               paste0("45_and_over_years:",ti$`45 and over years`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 248: GSE88824
{
  i=48
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("person_id:",ti$`person id`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 249:GSE86961
{
  i=49
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`tissue type`[r], # Sample type
                               ti$`subject status`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("race:",ti$race[r],";",
                                      "variant:",ti$variant[r],";",
                                      "multicentricity:",ti$multicentricity[r],";",
                                      "invasion:",ti$invasion[r],";",
                                      "nodal_status:",ti$nodal_status[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 250:GSE99994
{
  i=50
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("material:",ti$material[r],";",
                                      "tumor_subgroup:",ti$`tumor subgroup`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 251: GSE99996
{
  i=51
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("material:",ti$material[r],";",
                                      "tumor_subgroup:",ti$`tumor subgroup`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 252: GSE80762
{
  i=52
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("time_point:",ti$`time point`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 253:GSE74738
{
  i=53
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample tissue`[r], # Sample type
                               "NA", # Disease state
                               ti$`fetal sex/sex`[r], # Gender
                               ti$`age (years)`[r], # Age
                               "NA", # Anatomic Location
                               paste0("status_group:",ti$`status/group`[r],";",
                                      "age_info:years;",
                                      "fetal_gestational_age_weeks:",ti$`fetal gestational age (weeks)/trimester`[r],";",
                                      "plate_id:",ti$`450k plate id`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 254: GSE60274
{
  i=54
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("survival_status:",ti$survival_status[r],";",
                                      "geo_expression_data:",ti$`geo expression data`[r],";",
                                      "mgmt_status:",ti$mgmt_status[r],";",
                                      "growth_protocol:",ti$`growth protocol`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 255: GSE63409
{
  i=55
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               ti$`subject status`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("subject_id:",ti$`subject id`[r],";",
                                      "phenotype:",ti$phenotype[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 256: GSE73412
{
  i=56
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "blood", # Anatomic Location
                               paste0("origin:",ti$origin[r],";",
                                      "family_relationship:",ti$`family relationship`[r],";",
                                      "y_haplotype:",ti$`y haplotype`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 257: GSE89852
{
  i=57
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("hepatitis_virus:",ti$`hepatitis virus`[r],";",
                                      "age_by_decade:",ti$`age by decade`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 258: GSE65205
{
  i=58
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$asthma[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "nasal_epithelium", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 259: GSE94943
{
  i=59
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample type`[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 260:GSE104471
{
  i=60
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[t], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("sample_id:",ti$`sample id`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 261:GSE104472
{
  i=61
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("sample_id:",ti$`sample id`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 262: GSE104778
{
  i=62
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "cord_blood", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 263:GSE61450
{
  i=63
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$`subject status`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$`subject age`[r], # Age
                               "NA", # Anatomic Location
                               paste0("bmi:",ti$bmi[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 264: GSE85042
{
  i=64
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("age_of_mother:",ti$`age of the mother`[r],";",
                                      "smoking_cigarettes_during_pregnancy:",ti$`smoking cigarettes during pregnancy`[r],";",
                                      "gestational_stress:",ti$`gestational stress`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 265: GSE67444
{
  i=65
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("sample_blood_lead_level:",ti$`sample blood lead level`[r],";",
                                      "mothers_age_months:",ti$`mothers age (months)`[r],";",
                                      "gestational_age_months:",ti$`gestational age (months)`[r],";",
                                      "smoking:",ti$smoking[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 266:GSE89181
{
  i=66
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$histology[r], # Sample type
                               ti$`be status`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("stage:",ti$Stage[r],";",
                                 "run:",ti$run[r],";",
                                      "bmi:",ti$bmi[r],";",
                                      "smoke_cigarettes:",ti$smokecigarettes[r],";",
                                      "alcoholic_drinks_per_week:",ti$alcoholicrinksperweek[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 267: GSE98056
{
  i=67
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`tissue/cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               "NA", # Age
                               "blood", # Anatomic Location
                               paste0("status:",ti$status[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 268:GSE99511
{
  i=68
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample type`[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "cervico_vaginal_material", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 269:GSE103413
{
  i=69
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               ti$tissue[r], # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 270:GSE103911
{
  i=70
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("timepoint:",ti$timepoint[r],";",
                                      "future_cscc:",ti$`future cscc`[r],";",
                                      "current_cscc:",ti$`current cscc`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 271:GSE52025
{
  i=71
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "skin", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 272:GSE94785
{
  i=72
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$`age (years)`[r], # Age
                               "NA", # Anatomic Location
                               paste0("age_info:years;",
                                 "asbestos_exposed:",ti$`asbestos exposed`[r],";",
                                      "one_year_smoking_20_cigarettes_per_day:",ti$`one year of smoking 20 cigarettes/day)`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 273:GSE66210
{
  i=73
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("individual:",ti$individual[r],";",
                                      "pregnancy:",ti$pregnancy[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 274:GSE62219
{
  i=74
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$`age (months)`[r], # Age
                               "NA", # Anatomic Location
                               paste0("age_info:months")),# misc
                             nrow=1))
    message(r)
  }
}

# 275:GSE73895
{
  i=75
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "cancer", # Disease state
                               ti$gender[r], # Gender
                               ti$`subject age`[r], # Age
                               "glioblastoma", # Anatomic Location
                               paste0("alive:",ti$alive[r],";",
                                      "survival_time_months:",ti$`survival time (months)`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 276:GSE77269
{
  i=76
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$`subject status`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("patient_id:",ti$`patient id`[r],";",
                                      "tumor_size_cm:",ti$`tumor size (cm)`[r],";",
                                      "hbv_infection:",ti$`hbv infection`[r],";",
                                      "cirrhosis:",ti$cirrhosis[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 277:GSE72254
{
  i=77
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "breast_tumor", # Anatomic Location
                               paste0("cohort_id:",ti$`cohort id`[r],";",
                                      "grade:",ti$grade[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 278:GSE103769
{
  i=78
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("time:",ti$time[r],";",
                                      "outcome:",ti$outcome[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 279:GSE66313
{
  i=79
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$`subject age`[r], # Age
                               "NA", # Anatomic Location
                               paste0("future_development_of_breast_cancer:",ti$`future development of invasive breast cancer`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 280: GSE109042
{
  i=80
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$`disease status`[r], # Disease state
                               ti$Sex[r], # Gender
                               ti$age[r], # Age
                               "buccal_epithelium", # Anatomic Location
                               paste0("individual:",ti$individual[r],";",
                                      "self_reported_ethnicity:",ti$`self-reported ethnicity`[r],";",
                                      "disease_subgroup:",ti$`disease subgroup`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 281:GSE107352
{
  i=81
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 282: GSE61107
{
  i=82
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  dxvar <- ifelse(ti$`disease status (1=control, 2=scz patient)`==" 1","control", "SCZ patient")
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               dxvar[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 283: GSE89776
{
  i=83
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("sample_group:",ti$`sample  group`[r],";",
                                      "subject:",ti$subject[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 284:GSE76503
{
  i=84
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               ti$condition[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("intervention:",ti$intervention[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 285: GSE104287
{
  i=85
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "blood", # Anatomic Location
                               paste0("patient_id:",ti$`patient id`[r],";",
                                      "time:",ti$time[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 286:GSE102970
{
  i=86
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               "NA", # Age
                               "sperm", # Anatomic Location
                               paste0("tissue_subtype:",ti$`tissue subtype`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 287:GSE104812
{
  i=87
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$`age (y)`[r], # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 288:GSE98224
{
  i=88
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$diagnosis[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "placenta", # Anatomic Location
                               paste0("maternal_age:",ti$`maternal age`[r],";",
                                      "maternal_ethnicity:",ti$`maternal ethnicity`[r],";",
                                      "maternal_blood_type:",ti$`maternal blood type`[r],";",
                                      "fetal_weight_zscore:",ti$`fetal weight z-score`[r],";",
                                      "placenta_weight_zscore:",ti$`placenta weight z-score`[r],";",
                                      "mod:",ti$mod[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 289:GSE112812
{
  i=89
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("stability:",ti$stability[r],";",
                                      "group:",ti$group[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 290:GSE112873
{
  i=90
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("group:",ti$group[r],";",
                                      "stability:",ti$stability[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 291:GSE101641
{
  i=91
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$diagnosis[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "nasal_epithelium", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 292:GSE115797
{
  i=92
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample type`[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "skin", # Anatomic Location
                               paste0("batch:",ti$batch[r],";",
                                      "age_of_psoriasis_onset:",ti$`age of psoriasis onset`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 293:GSE85506
{
  i=93
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$group[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "blood", # Anatomic Location
                               paste0("inhibition_average:",ti$`inhibition (average values)`[r],";",
                                      "facilitation_average:",ti$`facilitation values (average)`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 294:GSE71719
{
  i=94
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("dna_type:",ti$`dna type`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 295: GSE56598
{
  i=95
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("tissue_type:",ti$`tissue type`[r],";",
                                      "subtype:",ti$subtype[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 296:GSE66552
{
  i=96
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "blood", # Anatomic Location
                               paste0("group:",ti$group[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 297:GSE104376
{
  i=97
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "cord_blood", # Anatomic Location
                               paste0("pregnancy_anxiety:",ti$`pregnancy anxiety`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 298:GSE119684
{
  i=98
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "blood", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 299:GSE116300
{
  i=99
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               ti$`case status`[r], # Disease state
                               ti$Sex[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("phenotype:",ti$phenotype[r],";",
                                      "mutation:",ti$mutation[r],";",
                                      "variant_classification:",ti$`variant classification`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 300:GSE104728
{
  i=100
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "medulloblastoma", # Anatomic Location
                               paste0("sample_type:",ti$`sample type`[r],";",
                                      "tumor_subgroup:",ti$`tumor subgroup`[r])),# misc
                             nrow=1))
    message(r)
  }
}


save(gat,file="geoanno_300gse.rda")
write.csv(gat,file="geoanno_300gse.csv")

#-------
# rep 7
#-------

length(unique(gat[,2])) # 300
tlsf <- tls[!names(tls) %in% unique(gat[,2])] # subset
length(tlsf) # 178
tlsf <- tlsf[sample(length(tlsf),length(tlsf))] # randomize order

# 301: GSE100940
{
  i=1
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 302:GSE117852
{
  i=2
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("genotype:",ti$genotype[r],";",
                                      "pannet_subtype:",ti$`pannet subtype`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 303: GSE75546
{
  i=3
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`tissue subtype`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "rectum", # Anatomic Location
                               paste0("node_metastasis_state:",ti$`node metastasis state`[r],";",
                                      "t:",ti$t[r],";",
                                      "n:",ti$n[r],";",
                                      "m:",ti$m[r],";",
                                      "stage:",ti$Stage[r],";",
                                      "differentiation:",ti$differentiation[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 304:GSE105798
{
  i=4
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 305: GSE74877
{
  i=5
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 306:GSE115920
{
  i=6
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$`age(years)`[r], # Age
                               "NA", # Anatomic Location
                               paste0("age_info:years;",
                                      "ethnicity:",ti$ethnicity[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 307:GSE86260
{
  i=7
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "prostate", # Anatomic Location
                               paste0("patient:",ti$patient[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 308:GSE79329
{
  i=8
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 309:GSE99652
{
  i=9
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "differentiation_status:",ti$`differentiation status`[r],";",
                                      "treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 310: GSE67419
{
  i=10
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`brain region`[r], # Sample type
                               ti$disease[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "brain", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 311: GSE65058
{
  i=11
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$status[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "liver", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 312:GSE89472
{
  i=12
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$`disease status`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$`age (years)`[r], # Age
                               "NA", # Anatomic Location
                               paste0("age_info:years;",
                                      "twin_pair:",ti$`twin pair`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 313:GSE62053
{
  i=13
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 314: GSE85845
{
  i=14
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               ti$`disease state`[r], # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("subject:",ti$subject[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 315:GSE107215
{
  i=15
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`developmental stage`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 316: GSE102504
{
  i=16
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$Sex[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 317:GSE108058
{
  i=17
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tag[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "sperm", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 318:GSE94962
{
  i=18
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("condition:",ti$condition[r],";",
                                      "genotype_variation:",ti$`genotype/variation`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 319:GSE79366
{
  i=19
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 320:GSE87154
{
  i=20
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell line`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("treatment:",ti$treatment[r],";",
                                      "control_vendor:",ti$`control vendor`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 321:GSE97483
{
  i=21
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`sample type`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 322:GSE63669
{
  i=22
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "cancer", # Disease state
                               ti$Sex[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("tumor_type:",ti$`tumor type`[r],";",
                                      "sample_id:",ti$`sample id`[r],";",
                                      "patient_id:",ti$`patient id`[r],";",
                                 "sample:",ti$sample[r],";",
                                      "histology:",ti$histology[r],";",
                                 "death:",ti$death[r],";",
                                 "subgroup:",ti$subgroup[r]
                                      )),# misc
                             nrow=1))
    message(r)
  }
}

# 323:GSE116754
{
  i=23
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("tissue:",ti$tissue[r],";",
                                      "developmental_stage:",ti$`developmental stage`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 324:GSE111942
{
  i=24
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "blood", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 325:GSE105067
{
  i=25
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "treatment:",ti$treatment[r],";",
                                      "day:",ti$day[r],";",
                                      "msssi:",ti$msssi[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 326:GSE54939
{
  i=26
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cdkn2a_status:",ti$`cdkn2a status`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 327: GSE78875
{
  i=27
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "treated_with:",ti$`treated with`[r],";",
                                      "subtype:",ti$subtype[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 328:GSE62178
{
  i=28
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 329: GSE90015
{
  i=29
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 330:GSE74167
{
  i=30
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`developmental stage`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 331:GSE69118
{
  i=31
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "resistant_phenotype:",ti$`resistant phenotype`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 332:GSE81224
{
  i=32
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$`tumor type`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$`patient age`[r], # Age
                               "NA", # Anatomic Location
                               paste0("tumor_histology:",ti$`tumor histology`[r],";",
                                      "stage:",ti$Stage[r],";",
                                      "grade:",ti$grade[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 333:GSE76372
{
  i=33
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 334:GSE80241
{
  i=34
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 335:GSE93208
{
  i=35
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell population`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               ti$tissue[r], # Anatomic Location
                               paste0("sample_type:",ti$`sample type`[r],";",
                                 "gestational_age_weeks:",ti$`gestational age (weeks)`[r],";",
                                      "date:",ti$date[r],";",
                                      "dna_quantity_ng:",ti$`dna quantity (ng)`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 336: GSE95486
{
  i=36
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "blood", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 337:GSE92911
{
  i=37
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               ti$`disease subtype`[r], # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 338:GSE65306
{
  i=38
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               ti$diagnosis[r], # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 339:GSE107226
{
  i=39
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               "NA", # Gender
                               ti$age[r], # Age
                               "lung", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 340:GSE92580
{
  i=40
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "brain", # Anatomic Location
                               paste0("dna_restoration:",ti$`dna restoration`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 341:GSE69634
{
  i=41
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "agent:",ti$agent[r],";",
                                      "dose:",ti$dose[r],";",
                                      "sample_type:",ti$`sample type`[r],";",
                                      "ip_antibody:",ti$`ip antibody`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 342:GSE87655
{
  i=42
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("stimulation:",ti$stimulation[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 343:GSE70739
{
  i=43
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 344:GSE117050
{
  i=44
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("patient:",ti$patient[r],";",
                                      "timepoint:",ti$timepoint[r],";",
                                      "future_cscc:",ti$`future cscc`[r],";",
                                      "current_cscc:",ti$`current cscc`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 345:GSE100563
{
  i=45
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell population`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               ti$`cell type`[r], # Anatomic Location
                               paste0("ethnic_group:",ti$`ethnic group`[r],";",
                                      "individual:",ti$individual[r],";",
                                      "infection:",ti$infection[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 346:GSE80377
{
  i=46
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "PBMC", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 347:GSE67477
{
  i=47
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "liver", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "agent:",ti$agent[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 348:GSE92843
{
  i=48
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 349:GSE89648
{
  i=49
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "cancer", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("tumor_stage:",ti$`tumor stage`[r],";",
                                      "passage_number_in_vitro:",ti$`passage number in vitro`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 350:GSE77206
{
  i=50
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "transfection:",ti$transfection[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 351: GSE85464
{
  i=51
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`tissue subtype`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "gastric", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 352:GSE101673
{
  i=52
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell or tissue type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("treated_with:",ti$`treated with`[r],";",
                                      "time_point:",ti$`time point`[r],";",
                                      "genotype_variation:",ti$`genotype/variation`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 353:GSE67350
{
  i=53
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "transgene:",ti$transgene[r],";",
                                      "treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 354:GSE115399
{
  i=54
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("gestational_age:",ti$`gestational age`[r],";",
                                      "exposure_duration:",ti$`exposure duration`[r],";",
                                      "exposed_to:",ti$`exposed to`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 355:GSE93963
{
  i=55
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "anatomic_site_of_tumor_origin:",ti$`anatomic site of tumour origin`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 356:GSE99650
{
  i=56
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "differentiation_status:",ti$`differentiation status`[r],";",
                                      "treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 357:GSE86833
{
  i=57
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 358:GSE86650
{
  i=58
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 359:GSE94956
{
  i=59
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("condition:",ti$condition[r],";",
                                      "genotype_variation:",ti$`genotype/variation`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 360:GSE75133
{
  i=60
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 361:GSE101864
{
  i=61
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell or tissue type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("genotype_variation:",ti$`genotype/variation`[r],";",
                                      "treated_with:",ti$`treated with`[r],";",
                                      "time_point:",ti$`time point`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 362:GSE66562
{
  i=62
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               ti$tissue[r], # Anatomic Location
                               paste0("donor_status:",ti$`donor status`[r],";",
                                      "donor_id:",ti$`donor id`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 363:GSE81015
{
  i=63
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "retinoblastoma", # Anatomic Location
                               paste0("cell_type:",ti$`cell type`[r],";",
                                      "treatment:",ti$treatment[r],";",
                                      "time_of_infection:",ti$`time of infection`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 364:GSE83842
{
  i=64
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$`age (yrs)`[r], # Age
                               "lung", # Anatomic Location
                               paste0("age_info:years;",
                                      "smoking_status:",ti$`smoking status`[r],";",
                                      "brinkman_index:",ti$`brinkman index`[r],";",
                                      "tnm_factor:",ti$`tnm factor`[r],";",
                                      "stage:",ti$Stage[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 365:GSE117853
{
  i=65
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("genotype:",ti$genotype[r],";",
                                      "pannet_subtype:",ti$`pannet subtype`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 366:GSE52980
{
  i=66
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               ti$`body site`[r], # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 367:GSE95761
{
  i=67
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 368:GSE106600
{
  i=68
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("disease_state:",ti$`disease state`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 369:GSE92469
{
  i=69
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$age[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 370:GSE74214
{
  i=70
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               ti$`subject age`[r], # Age
                               "breast", # Anatomic Location
                               paste0("subject_bmi:",ti$`subject bmi`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 371:GSE86258
{
  i=71
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "prostate", # Anatomic Location
                               paste0("patient:",ti$patient[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 372:GSE89925
{
  i=72
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 373:GSE95488
{
  i=73
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "blood", # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "blood", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 374:GSE110607
{
  i=74
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$status[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("zygosity:",ti$zygosity[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 375:GSE107039
{
  i=75
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "liver", # Sample type
                               ti$disease[r], # Disease state
                               ti$Sex[r], # Gender
                               ti$age[r], # Age
                               "liver", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 376:GSE71525
{
  i=76
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample type`[r], # Sample type
                               "NA", # Disease state
                               ti$tag[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("sample_source_id:",ti$`sample source id`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 377:GSE107039
{
  i=77
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 378:GSE66459
{
  i=78
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "blood", # Anatomic Location
                               paste0("gestational_age_days:",ti$`gestational_age (days)`[r],";",
                                      "initiation_of_labor:",ti$`initiation of labor`[r],";",
                                      "")),# misc
                             nrow=1))
    message(r)
  }
}

# 379:GSE79064
{
  i=79
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               ti$`disease status`[r], # Disease state
                               ti$Sex[r], # Gender
                               "NA", # Age
                               ti$tissue[r], # Anatomic Location
                               paste0("age_of_disease_onset:",ti$`age of disease onset`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 380:GSE89401
{
  i=80
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 381: GSE103425
{
  i=81
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "vector:",ti$vector[r],";",
                                      "atra:",ti$atra[r],";",
                                      "decitabine:",ti$decitabine[r],";",
                                      "replicate:",ti$replicate[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 382:GSE87056
{
  i=82
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`tissue region`[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "liver", # Anatomic Location
                               paste0("individual:",ti$individual[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 383:GSE63267
{
  i=83
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 384:GSE105288
{
  i=84
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 385:GSE86829
{
  i=85
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 386:GSE97853
{
  i=86
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 387:GSE92462
{
  i=87
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$age[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 388:GSE116992
{
  i=88
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "blood", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 389:GSE75405
{
  i=89
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$`subject status`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$`age (yrs)`[r], # Age
                               "peripheral blood", # Anatomic Location
                               paste0("age_info:years;",
                                      "subject_id:",ti$`subject id`[r],";",
                                      "sorted_population:",ti$`sorted population`[r],";",
                                      "sorted_population_markers:",ti$`sorted population markers`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 390:GSE107038
{
  i=90
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "liver", # Sample type
                               ti$disease[r], # Disease state
                               ti$Sex[r], # Gender
                               ti$age[r], # Age
                               "liver", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 391: GSE79648
{
  i=91
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 392:GSE98813
{
  i=92
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("generations:",ti$generations[r],";",
                                      "treatment:",ti$treatment[r],";",
                                      "cetuximab_resistant_clone:",ti$`cetuximab resistant clone`[r],";",
                                      "batch:",ti$batch[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 393:GSE73950
{
  i=93
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$`menstrual phase`[r], # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("patient_id:",ti$`patient id`[r],";",
                                      "menstrual_phase:",ti$`menstrual phase`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 394:GSE91071
{
  i=94
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 395:GSE65057
{
  i=95
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$status[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "liver", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 396:GSE95058
{
  i=96
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("donor:",ti$donor[r],";",
                                      "passage:",ti$passage[r],";",
                                      "differentiation_condition:",ti$`differentiation condition`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 397:GSE82234
{
  i=97
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "umbilical_vein_endothelium", # Anatomic Location
                               paste0("donor:",ti$donor[r],";",
                                      "passage:",ti$passage[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 398:GSE106099
{
  i=98
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("fetal_sex:",ti$`fetal sex`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 399:GSE96864
{
  i=99
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell type`[r],";",
                                      "transfected_with:",ti$`transfected with`[r],";",
                                      "time_point:",ti$`time point`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 400: GSE61461
{
  i=100
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# study   
{
  i=
  ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

save(gat,file="geoanno_399gse.rda")
write.csv(gat,file="geoanno_399gse.csv")

#-------
# rep 8
#-------
load("geo_gse-atables_list.rda")
#sv <- unlist(lapply(tgse.list,nrow))
#tls <- tgse.list[rev(order(sv))]
tls <- tgse.list

length(unique(gat[,2])) # 399
tlsf <- tls[!names(tls) %in% unique(gat[,2])] # subset
length(tlsf) # 79
tlsf <- tlsf[sample(length(tlsf),length(tlsf))] # randomize order

# 401: GSE79100
{
  i=1
  ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               ti$age[r], # Age
                               "kidney", # Anatomic Location
                               paste0("ethnicity:",ti$ethnicity[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 402: GSE65078
{
  i=2
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 403: GSE108562
{
  i=3
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               ti$tissue[r], # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 404: GSE81846
{
  i=4
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 405: GSE74233
{
  i=5
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "stomach", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 406: GSE108564
{
  i=6
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               ti$tissue[r], # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 407: GSE85647
{
  i=7
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$`clinical phenotype`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$`age (years)`[r], # Age
                               "NA", # Anatomic Location
                               paste0("age_info:years")),# misc
                             nrow=1))
    message(r)
  }
}

# 408: GSE120307
{
  i=8
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$diagnosis[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "peripheral_blood", # Anatomic Location
                               paste0("pair_number:",ti$`pair number`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 409: GSE77353 
{
  i=9
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("histone_h3_mutations:",ti$`histone h3 mutations`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 410:GSE100249
{
  i=10
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 411: GSE115783
{
  i=11
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "pituitary_gland", # Anatomic Location
                               paste0("tumor_type:",ti$`tumor type`[r],";",
                                      "invasiveness:",ti$invasiveness[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 412: GSE75196
{
  i=12
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$disease[r], # Disease state
                               ti$Sex[r], # Gender
                               "NA", # Age
                               "placenta", # Anatomic Location
                               paste0("gestation_week:",ti$`gestation (wk)`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 413:GSE89649
{
  i=13
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "cancer", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("passage_number_in_vitro:",ti$`passage number in vitro`[r],";",
                                      "tumor_stage:",ti$`tumor stage`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 414:GSE107511
{
  i=14
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample type`[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               ti$tissue[r], # Anatomic Location
                               paste0("tumor_subgroup:",ti$`tumor subgroup`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 415:GSE101658
{
  i=15
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$diagnosis[r], # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "hippocampus", # Anatomic Location
                               paste0("myelination_status:",ti$`myelination status`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 416:GSE59524
{
  i=16
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 417: GSE95036
{
  i=17
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               ti$`anatomical site`[r], # Anatomic Location
                               paste0("tissue_storage:",ti$`tissue storage`[r],";",
                                      "center:",ti$center[r],";",
                                      "hpv_status:",ti$`hpv status`[r],";",
                                      "organ_group:",ti$`organ group`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 418: GSE81006
{
  i=18
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "genotype:",ti$genotype[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 419:GSE94063
{
  i=19
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "sample_group:",ti$`sample groupd`[r],";",
                                      "treatment:",ti$treatment[r],";",
                                      "time:",ti$time[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 420:GSE64096
{
  i=20
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("gradient_layer:",ti$`gradient layer`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 421:GSE67097
{
  i=21
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               ti$`body site`[r], # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 422:GSE75406
{
  i=22
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$`subject status`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$`age (yrs)`[r], # Age
                               ti$tissue[r], # Anatomic Location
                               paste0("subject_id:",ti$`subject id`[r],";",
                                      "age_info:years;",
                                      "sorted_population:",ti$`sorted population`[r],";",
                                      "sorted_population_markers:",ti$`sorted population markers`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 423:GSE86829
{
  i=23
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 424:GSE118570
{
  i=24
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 425: GSE90871
{
  i=25
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "prefrontal_cortex", # Anatomic Location
                               paste0("maternal_smoking:",ti$`maternal smoking`[r],";",
                                      "stage:",ti$Stage[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 426:GSE72338
{
  i=26
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("subject_id:",ti$`subject id`[r],";",
                                      "study_group:",ti$`study group`[r],";",
                                      "matched_pair_id:",ti$`matched pair id`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 427:GSE82084
{
  i=27
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               "NA", # Age
                               "cord_blood", # Anatomic Location
                               paste0("gestational_age:",ti$`gestational age`[r],";",
                                      "individual:",ti$individual[r],";",
                                      "tissue:",ti$tissue[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 428:GSE74609
{
  i=28
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               gsub("[a-z|A-Z]| ","",ti$age[r]), # Age
                               "bone_marrow", # Anatomic Location
                               paste0("group:",ti$group[r],";",
                                      "race:",ti$race[r],";",
                                      "time_point:",ti$`time point`[r],";",
                                      "age_info:years")),# misc
                             nrow=1))
    message(r)
  }
}

# 429: GSE95061
{
  i=29
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("donor:",ti$donor[r],";",
                                      "passage:",ti$passage[r],";",
                                      "differentiation_condition:",ti$`differentiation condition`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 430:GSE97484
{
  i=30
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`sample type`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 431: GSE60655
{
  i=31
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "thigh_muscle", # Anatomic Location
                               paste0("batch:",ti$batch[r],";",
                                      "group:",ti$group[r],";",
                                      "subject:",ti$subject[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 432:GSE65079
{
  i=32
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 433:GSE107737
{
  i=33
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "whole_blood", # Sample type
                               ti$diagnosis[r], # Disease state
                               ti$Sex[r], # Gender
                               gsub(" |[a-z|A-Z]","",ti$age[r]), # Age
                               ti$tissue[r], # Anatomic Location
                               paste0("age_info:years")),# misc
                             nrow=1))
    message(r)
  }
}

# 434: GSE87582
{
  i=34
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample type`[r], # Sample type
                               ti$group[r], # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("cell_type:",ti$`cell type`[r],";",
                                      "education:",ti$education[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 435: GSE67351
{
  i=35
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "transgene:",ti$transgene[r],";",
                                      "treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 436: GSE73948
{
  i=36
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "uterus", # Anatomic Location
                               paste0("patient_id:",ti$`patient id`[r],";",
                                      "menstrual_phase:",ti$`menstrual phase`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 437: GSE94282
{
  i=37
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "sample_group:",ti$`sample groupd`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 438: GSE99716
{
  i=38
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_source:",ti$`cell source`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 439:GSE80794
{
  i=39
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "breast_cancer", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 440:GSE62670
{
  i=40
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 441: GSE120062
{
  i=41
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "placenta", # Anatomic Location
                               paste0("weight:",ti$weight[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 442: GSE77965
{
  i=42
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`tissue type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               ti$tissue[r], # Anatomic Location
                               paste0("individual:",ti$individual[r],";",
                                      "sample_group1:",ti$samplegroup[r],";",
                                      "sample_group2:",ti$`sample group`[r],";",
                                      "xenograft:",ti$xenograft[r],";",
                                      "cell_line:",ti$`cell line`[r],";",
                                      "treatment:",ti$treatment[r],";",
                                      "")),# misc
                             nrow=1))
    message(r)
  }
}

# 443:GSE81308
{
  i=43
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 444:GSE113779
{
  i=44
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "esophagus", # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "esophagus", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 445:GSE75550
{
  i=45
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`tissue subtype`[r], # Sample type
                               "cancer", # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "rectum", # Anatomic Location
                               paste0("subject_id:",ti$`subject id`[r],";",
                                      "stage:",ti$Stage[r],";",
                                      "node_metastasis_state:",ti$`node metastasis state`[r],";",
                                      "differentiation:",ti$differentiation[r],";",
                                      "t:",ti$t[r],";",
                                      "n:",ti$n[r],";",
                                      "m:",ti$m[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 446: GSE105093
{
  i=46
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample type`[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("reprogramming_type:",ti$`reprogramming type`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 447:GSE89400
{
  i=47
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "glioma", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "cell_type:",ti$`cell type`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 448: GSE73626
{
  i=48
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample type`[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$`age (yrs)`[r], # Age
                               ti$joint[r], # Anatomic Location
                               paste0("age_info:years;",
                                      "body_mass_index:",ti$`body mass index`[r],";",
                                      "subject_status:",ti$`subject status`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 449:GSE85942
{
  i=49
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("genotype:",ti$genotype[r],";",
                                      "dox:",ti$dox[r],";",
                                      "passage:",ti$passage[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 450:GSE101443
{
  i=50
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`tumor status`[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("patient_id:",ti$`patient id`[r],";",
                                      "tet1_qpcr_level:",ti$`tet1 qpcr level`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 451:GSE78732
{
  i=51
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("morphology:",ti$morphology[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 452:GSE77036
{
  i=52
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 453:GSE80017
{
  i=53
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$`disease status`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "prefrontal_cortex", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 454:GSE113061
{
  i=54
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 455:GSE112012
{
  i=55
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 456:GSE87797
{
  i=56
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("donor:",ti$donor[r],";",
                                      "passage:",ti$passage[r],";",
                                      "group:",ti$group[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 457:GSE78279
{
  i=57
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("transformation_stage:",ti$`transformation stage`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 458:GSE100653
{
  i=58
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               ti$`subject age`[r], # Age
                               "breast", # Anatomic Location
                               paste0("subject_bmi:",ti$`subject bmi`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 459:GSE101445
{
  i=59
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`tumor status`[r], # Sample type
                               "cancer", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("patient_id:",ti$`patient id`[r],";",
                                      "tet1_qpcr_level:",ti$`tet1 qpcr level`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 460:GSE87798
{
  i=60
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("donor:",ti$donor[r],";",
                                      "group:",ti$group[r],";",
                                      "passage:",ti$passage[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 461:GSE73949
{
  i=61
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               ti$condition[r], # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "uterus", # Anatomic Location
                               paste0("patient_id:",ti$`patient id`[r],";",
                                      "menstrual_phase:",ti$`menstrual phase`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 462:GSE81025
{
  i=62
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 463:GSE73745
{
  i=63
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample type`[r], # Sample type
                               ti$`disease state`[r], # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 464: GSE98938
{
  i=64
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$tissue[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("gestational_stage:",ti$`gestational stage`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 465: GSE100503
{
  i=65
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("matched_surgical_specimen_accession_gse66313:",ti$`matched surgical specimen accession (gse66313)`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 466:GSE76709
{
  i=66
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 467:GSE94350
{
  i=67
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("treatment:",ti$treatment[r],";",
                                      "time:",ti$time[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 468:GSE77136
{
  i=68
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell line`[r], # Sample type
                               "NA", # Disease state
                               ti$Sex[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("race:",ti$race[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 469:GSE69852
{
  i=69
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               ti$age[r], # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 470:GSE90060
{
  i=70
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$type[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               ti$age[r], # Age
                               "uterus", # Anatomic Location
                               paste0("time_point:",ti$`time point`[r],";",
                                      "menstrual_phase:",ti$`menstrual phase`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 471:GSE74797
{
  i=71
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 472:GSE91069
{
  i=72
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 473:GSE83261
{
  i=73
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "NA", # Disease state
                               ti$gender[r], # Gender
                               "NA", # Age
                               "skin", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 474:GSE85649
{
  i=74
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               ti$`clinical phenotype`[r], # Disease state
                               ti$gender[r], # Gender
                               ti$`age (years)`[r], # Age
                               "NA", # Anatomic Location
                               paste0("age_info:years")),# misc
                             nrow=1))
    message(r)
  }
}

# 475:GSE98815
{
  i=75
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample type`[r], # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("batch:",ti$batch[r],";",
                                      "generations:",ti$generations[r],";",
                                      "treatment:",ti$treatment[r],";",
                                      "cetuximab_resistant_clone:",ti$`cetuximab resistant clone`[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 476:GSE81790
{
  i=76
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("treatment:",ti$treatment[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 477:GSE103427
{
  i=77
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "cell_line", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("cell_line:",ti$`cell line`[r],";",
                                      "vector:",ti$vector[r],";",
                                      "atra:",ti$atra[r],";",
                                      "decitabine:",ti$decitabine[r],";",
                                      "replicate:",ti$replicate[r])),# misc
                             nrow=1))
    message(r)
  }
}

# 478:GSE68777
{
  i=78
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`sample type`[r], # Sample type
                               ti$diagnosis[r], # Disease state
                               ti$Sex[r], # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}

# 479:GSE63670
{
  i=79
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               ti$`cell type`[r], # Sample type
                               "cancer", # Disease state
                               ti$Sex[r], # Gender
                               ti$age[r], # Age
                               "brain", # Anatomic Location
                               paste0("sample:",ti$sample[r],";",
                                      "sample_id:",ti$`sample id`[r],";",
                                      "patient_id:",ti$`patient id`[r],";",
                                      "tumor_type:",ti$`tumor type`[r],";",
                                      "histology:",ti$histology[r],";",
                                      "death:",ti$death[r],";",
                                      "subgroup:",ti$subgroup[r])),# misc
                             nrow=1))
    message(r)
  }
}

# study   
{
  i=
    ti <- as.data.frame(tlsf[[i]], stringsAsFactors = F)
  head(ti)
  for(r in 1:nrow(ti)){
    gat <- rbind(gat, matrix(c(ti$gsm[r], # GSM ID
                               ti$gse[r], # GSE ID
                               "NA", # Sample type
                               "NA", # Disease state
                               "NA", # Gender
                               "NA", # Age
                               "NA", # Anatomic Location
                               paste0("NA")),# misc
                             nrow=1))
    message(r)
  }
}


save(gat,file="geoanno_478gse.rda")
write.csv(gat,file="geoanno_478gse.csv")