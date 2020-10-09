#!/usr/bin/env/R

# Describes postprocessing from coerced GEO GSM metadata.
# Notes on rules for regex matching:
# 1. var states having phrases separated by spaces, to be separated by underscores
# 2. variable with more than one state have states separated by semicolons
# 3. all var states tolower (lowercase only)
# 4. auto-populate redundant variable states
# 5. add negative match check for disease status (e.g. excludes "non-cancer" from "cancer" search)

library(data.table)

load("md-preprocess.rda") # coerced, partially annotated metadata
load("mdmap-gsm_35k.rda") # MetaSRA-pipeline, mapped and predicted labels
ccf = fread('ccformat.txt', sep = ' ', header = T) # formatted Cellosaurus records
load("prepmd.rda") # storage procedure annotations

md = md.preprocess
mdpost = md[,c(1, 2, 3)]
nfn = "md-postprocess.rda" # new file name
mdpost$sampletype = mdpost$tissue = mdpost$disease = "NA"
mdpost$arrayid_full = paste0(md$array_id, "_", md$sentrix_id)
mdpost$basename = md$basename

#-----------------
# helper functions
#-----------------
get_pstr = function(v){
  # get automatic regex patterns
  # does progressive capitalization on values separated by spaces
  # for each value in v, appends flanking '.*' (matches any char)
  # for each value in v, appends "|" OR conditional separator
  # nfilt: Boolean, adds negative statement filter.
  rs = ""
  
  for(ci in 1:length(v)){
    c = v[ci]
    if(ci == 1){
      rs = paste0(".*", c, ".*")
    } else{
      rs = paste(c(rs, paste0(".*", c, ".*")), collapse = "|")
    }
    uv = unlist(strsplit(c, " ")) # num space-sep units
    # for each unit use lower- and uppercase
    uvstr = c(paste(uv, collapse = " "))
    uvl = list(uv)
    for(i in 1:length(uv)){
      uvi = c()
      for(ui in 1:length(uv)){
        chari = uv[ui]
        if(ui <= i){
          if(nchar(chari)>1){
            ssi = paste0(toupper(substr(chari, 1, 1)),
                         substr(chari, 2, nchar(chari)))
          } else{
            ssi = paste0(toupper(substr(chari, 1, 1)))
          }
        }
        else{
          ssi = chari
        }
        uvi = c(uvi, ssi)
      }
      uvl[[length(uvl)+1]] = uvi
    }
    # append to new str
    for(si in 1:length(uvl)){
      s = uvl[[si]]
      if(length(uv) > 1){
        if(!si==1){
          # space sep
          rs = paste(c(rs, paste0(".*", paste(s, collapse = " "), ".*")), collapse = "|")
        }
        # underline sep
        rs = paste(c(rs, paste0(".*", paste(s, collapse = "_"), ".*")), collapse = "|")
        # dash sep
        rs = paste(c(rs, paste0(".*", paste(s, collapse = "-"), ".*")), collapse = "|")
      } else{
        if(!si==1){
          rs = paste(c(rs, paste0(".*", s, ".*")), collapse = "|")
        }
      }
    }
  }
  return(rs)
}

get_pstr_neg = function(pstr){
  # pstr: output of get_pstr
  # returns patterns for exclusion/negative lookup
  pstrg = gsub("\\.\\*", "_", pstr); pstrg = gsub("\\|", "", pstrg)
  uv = unlist(strsplit(pstrg, "_")); uv = uv[!uv==""]
  for(ui in 1:length(uv)){
    s = uv[ui]
    if(ui == 1){
      ns = paste(paste0(".*non-", s, ".*"), paste0(".*Non-", s, ".*"), paste0(".*non ", s, ".*"),
                 paste0(".*Non ", s, ".*"), paste0(".*not ", s, ".*"), paste0(".*Not ", s, ".*"),
                 paste0(".*", s, "-free.*"), paste0(".*", s, "-Free.*"), paste0(".*", s, " free.*"),
                 paste0(".*", s, " Free.*"), sep = "|")
    } else{
      ns = paste(ns, paste0(".*non-", s, ".*"), paste0(".*Non-", s, ".*"), paste0(".*non ", s, ".*"),
                 paste0(".*Non ", s, ".*"), paste0(".*not ", s, ".*"), paste0(".*Not ", s, ".*"),
                 paste0(".*", s, "-free.*"), paste0(".*", s, "-Free.*"), paste0(".*", s, " free.*"),
                 paste0(".*", s, " Free.*"), sep = "|")
    }
  }
  return(ns)
}

get_filt = function(v, filtrel = "|",
                    varl = c("gsm_title", "sample_type", "disease_state", "anatomic_location", "misc"),
                    nfilt = FALSE, ntfilt = "", ptfilt = "", m = md){
  # Returns vector of conditionals corresponding to pattern match across vars in m
  # v: Reg. ex. output of get_pstr
  # nfilt: String output from 'get_pstr_neg', specifying negative lookup/exclusions
  # tfilt: Character vector, excludes records on term match from 'get_pstr'
  
  if(!filtrel %in% c("|", "&")){
    message("Please provide a valid filter relation symbol.")
    return(NULL)
  }
  # positive match filter
  if(ptfilt == ""){
    filtl = grepl(v, m[,colnames(m)==varl[1]])
  } else{
    filtl = grepl(get_pst(ptfilt), m[,colnames(m)==varl[1]])
    filtl = grepl(v, m[,colnames(m)==varl[1]])
  }
  # negative match filter
  if(nfilt){
    message("Using negative lookup/exclusion filter...")
    nfiltv = get_pstr_neg(v)
    filtl = filtl & !grepl(nfiltv, m[,colnames(m)==varl[1]])
  }
  # term filter
  if(!ntfilt == ""){
    message("Using term lookup filter...")
    filtl = filtl & !grepl(get_pstr(ntfilt), m[,colnames(m)==varl[1]])
  }
  # proceed if additional vars specified
  if(length(varl)>1){
    for(vi in varl[2:length(varl)]){
      if(filtrel == "|"){
        filtl = filtl | grepl(v, m[,colnames(m)==vi])
        if(nfilt){
          filtl = filtl & !grepl(nfiltv, m[,colnames(m)==vi])
        }
      }
      if(filtrel == "&"){
        filtl = filtl | grepl(v, m[,colnames(m)==vi])
        if(nfilt){
          filtl = filtl & !grepl(nfiltv, m[,colnames(m)==vi])
        }
      }
    }
  }
  return(filtl)
}

appendvar = function(var, val, filtv, m = mdpost){
  # Returns the var in m with the appended value
  # Replaces NA terms, append new terms, does not append repeated terms
  varr = m[, colnames(m) == var]
  # get composite filter
  filti = !grepl(paste0("(^|;)", val, "(;|$)"), varr); compfilt = filti & filtv
  # assess filter results
  if(length(compfilt[compfilt]) == 0){
    message("No unique values to append. Returning var unchanged.")
    return(varr)
  } else{
    varr[compfilt] = ifelse(varr[compfilt] == "NA", val,
                            paste(varr[compfilt], val, sep = ";")
    )
    message("Appended n = ", length(varr[compfilt]), " values")
    return(varr)
  }
}

#---------------
# disease status
#---------------
# Note: term matching also excludes negative matches with 'nfilt = T'
# e.g. non-cancer, not cancer, Not cancer, cancer-free, etc.
whichvar = c("gsm_title", "tissuevar", "disease_state", "anatomic_location") # colnames of variables to search in md
{
  # general terms
  {
    # chronic
    ssv = c("chronic"); sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "chronic", dfilt)
    # acute
    ssv = c("acute"); sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "acute", dfilt)
    # syndrome
    ssv = c("syndrome"); sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "syndrome", dfilt)
    # disorder
    ssv = c("disorder"); sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "disorder", dfilt)
    # inflammation
    ssv = c("inflam"); sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "inflammation", dfilt)
    # disorder
    ssv = c("disorder"); sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "disorder", dfilt)
  }

  # cancer, major terms
  ssv = c("cancer")
  sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
  mdpost$disease = appendvar("disease", "cancer", dfilt)

  # study groups, inc. control, case, healthy, etc.
  {
    # normal
    ssv = c("case"); sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "case", dfilt)
    # normal
    ssv = c("normal"); sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "normal", dfilt)
    # healthy
    ssv = c("healthy"); sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "healthy", dfilt)
    # control
    ssv = c("control"); sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "control", dfilt)
  }

  # Psychiatric and neurodegenerative
  {
    # psychiatric disorder
    ssv = c("alzheimer's", "alzheimers", "dementia", "anxi", "depression", "attention deficit", "ADHD")
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "alzheimers", dfilt)
    # alzheimer's
    ssv = c("alzheimer's", "alzheimers"); sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "alzheimers", dfilt)
    # dementia
    ssv = c("dementia"); sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "dementia", dfilt)
    # anxiety
    ssv = c("anxi"); sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "anxiety", dfilt)
    # depression
    ssv = c("depression"); sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "depression", dfilt)
    # attention deficit
    ssv = c("attention deficit", "ADHD"); sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "attention_deficit", dfilt)
  }

  # Arthritis, inc. fibromyalgia, gout, etc.
  {
    # arthritis
    ssv = c("arthritis", "rheumatoid", "psoriatic", "fibromyalgia", "gout");
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "arthritis", dfilt)
    # rheumatoid arthritis
    ssv = c("rheumatoid arthritis"); sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "rhematoid_arthritis", dfilt)

    # osteoarthritis
    ssv = c("osteoarthritis")
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "osteoarthritis", dfilt)

    # psoriatic arthritis
    ssv = c("psoriatic")
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "psoriatic_arthritis", dfilt)

    # fibromyalgia
    ssv = c("fibromyalgia")
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "fibromyalgia", dfilt)

    # gout
    ssv = c("gout")
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "gout", dfilt)
  }

  # Chromosomal abnormalities
  {
    ssv = c("trisomy", "monosomy", "triploidy",
            "chromosomal duplication", "chromosomal deletion",
            "down syndrome", "cri du chat",
            "balanced translocation", "unbalanced translocation",
            "pallister killian", "ring chromosome",
            "deletion syndrome", "klinefelter", "XXY", "turner", "mosaicism", "XXY")
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "chromosomal_abnormality", dfilt)
    # Down syndrome
    ssv = c("down syndrome", "trisomy 21", "trisomy twentyone", "trisomy twenty one")
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "down_syndrome", dfilt)
    # Klinefelter syndrome
    ssv = c("klinefelter", "XXY")
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "klinefelter_syndrome", dfilt)
  }

  # Genetic disorder
  {
    ssv = c("fragile x", "cystic fibrosis", "duane", "polycystic", "chrons",
            "hemophelia", "haemophelia", "hemochromatosis", "huntington's", "huntingtons",
            "thalassemia", "tay sachs", "tay sach", "parkinson's", "parkinsons",
            "sickle cell", "marfan")
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "genetic_disorder", dfilt)
    # Fragile X
    ssv = c("fragile x")
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "fragile_x", dfilt)
    # Cystic Fibrosis
    ssv = c("cystic fibrosis")
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "cystic_fibrosis", dfilt)
    # Thalassemia
    ssv = c("thalassemia")
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "thalassemia", dfilt)
    # Tay Sachs
    ssv = c("tay sachs", "tay sach")
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "tay_sachs", dfilt)
    # Parkinson's
    ssv = c("parkinson's", "parkinsons")
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "parkinsons_disease", dfilt)
    # Huntington's
    ssv = c("huntington's", "huntingtons")
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "huntingtons_disease", dfilt)
    # Sickle Cell Anemia
    ssv = c("sickle cell")
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "sickle_cell_anemia", dfilt)
    # Marfan
    ssv = c("marfan")
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "marfan_syndrome", dfilt)
    # Hemophelia
    ssv = c("hemophelia", "haemophilia")
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "hemophelia", dfilt)
  }

  # Other common conditions
  {
    # Anemia
    ssv = c("anemi")
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "anemia", dfilt)
    # Atrophy
    ssv = c("atroph")
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "atrophy", dfilt)
    # Scoliosis
    ssv = c("scoliosis")
    sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "scoliosis", dfilt)
    # Obesity
    ssv = c("obese"); sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "obese", dfilt)
    # Asthma
    ssv = c("asthma"); sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, nfilt = T, varl = whichvar)
    mdpost$disease = appendvar("disease", "asthma", dfilt)
  }

}

save(mdpost, file = nfn)

#------------------------------------
# tissue and disease, cancer subtypes
#------------------------------------
# Note: borrows heavily from types studied in TCGA
# Note: sped up by excluding 'nfilt = T' (less likely issue for specific-subtype details)
whichvar = c("gsm_title", "tissuevar", "disease_state", "anatomic_location")
{
  # Cancer tissue types
  {
    # tumor
    ssv = c("tumor", "tumour")
    mdpost$tissue = appendvar("tissue", "tumor", get_filt(get_pstr(ssv), nfilt = T, varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))
    # metastasis
    ssv = c("metasta")
    mdpost$tissue = appendvar("tissue", "metastasis", get_filt(get_pstr(ssv), nfilt = T, varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))
    # carcinoma
    ssv = c("carcinoma")
    mdpost$tissue = appendvar("tissue", "matched", get_filt(get_pstr(ssv), nfilt = T, varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))
    # sarcoma
    ssv = c("sarcoma")
    mdpost$tissue = appendvar("tissue", "matched", get_filt(get_pstr(ssv), nfilt = T, varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))
    # neoplasm
    ssv = c("neoplas")
    mdpost$tissue = appendvar("tissue", "matched", get_filt(get_pstr(ssv), nfilt = T, varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))
    # adenocarcinoma
    ssv = c("adenocarcinoma")
    mdpost$tissue = appendvar("tissue", "matched", get_filt(get_pstr(ssv), nfilt = T, varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))
  }

  # Cancer, disease subtypes
  {
    ssv = c("mesothelioma")
    mdpost$disease = appendvar("disease", "mesothelioma", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("melanoma", "skin cancer")
    mdpost$disease = appendvar("disease", "skin_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("glioblastoma", "glioma", "astrocytoma", "brain cancer")
    mdpost$disease = appendvar("disease", "brain_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("breast lobular carcinoma", "breast ductal carcinoma", "breast cancer")
    mdpost$disease = appendvar("disease", "breast_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("colorectal adeno", "colon cancer", "colorectal cancer", "rectal cancer")
    mdpost$disease = appendvar("disease", "colorectal_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("stomach adeno", "stomach cancer", "gastric cancer", "gastric adeno")
    mdpost$disease = appendvar("disease", "stomach_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("esophageal carcinoma", "esophageal adeno", "esophageal squamous cell carcinoma","oesophageal carcinoma", "oesophageal adeno", "oesophageal squamous cell carcinoma", "esophageal cancer", "oesophageal cancer")
    mdpost$disease = appendvar("disease", "esophageal_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("paraganglioma", "nerve cancer", "nerve cell cancer")
    mdpost$disease = appendvar("disease", "nerve_cell_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("cholangiocarcinoma", "bile duct cancer")
    mdpost$disease = appendvar("disease", "bile_duct_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("ovarian serous carcinoma", "ovarian epithelial cancer", "ovarian cancer")
    mdpost$disease = appendvar("disease", "ovarian_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("uterine carcinosarcoma", "uterine corpus endometrial carcinoma", "endometrial carcinoma", "uterine cancer")
    mdpost$disease = appendvar("disease", "uterine_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("cervical squamous cell carcinoma", "cervical squamous cell adenocarcinoma", "cervical cancer")
    mdpost$disease = appendvar("disease", "cervical_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("head and neck squamous cell carcinoma", "head and neck cancer")
    mdpost$disease = appendvar("disease", "head_and_neck_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("thyroid carcinoma", "thyroid cancer")
    mdpost$disease = appendvar("disease", "thyroid_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("lung adenocarcinoma", "lung squamous cell carcinoma", "lung cancer")
    mdpost$disease = appendvar("disease", "lung_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("clear cell renal cell carcinoma", "chromophobe renal cell carcinoma", "renal cancer", "kidney papillary carcinoma", "kidney cancer")
    mdpost$disease = appendvar("disease", "kidney_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("invasive urothelial bladder cancer", "bladder cancer")
    mdpost$disease = appendvar("disease", "bladder_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("prostate adenocarcinoma", "prostate cancer")
    mdpost$disease = appendvar("disease", "prostate_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("adrenocortical carcinoma", "pheochromocytoma", "adrenal cancer", "adrenal gland cancer")
    mdpost$disease = appendvar("disease", "adrenal_gland_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("liver hepatocellular carcinoma", "hepatoblastoma", "cholangiocarcinoma", "liver angiosarcoma", "liver cancer")
    mdpost$disease = appendvar("disease", "liver_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("pancreatic ductal adenocarcinoma", "pancreatic cancer")
    mdpost$disease = appendvar("disease", "pancreatic_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("uveal melanoma", "uveal lymphoma", "intraocular cancer", "retinoblastoma", "retinal cancer")
    mdpost$disease = appendvar("disease", "eye_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("thymoma", "thymus cancer", "thymus gland cancer", "thymic cancer")
    mdpost$disease = appendvar("disease", "thymus_gland_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))

    ssv = c("testicular germ cell cancer", "testicular cancer")
    mdpost$disease = appendvar("disease", "testicular_cancer", get_filt(get_pstr(ssv), varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), varl = whichvar))
  }

  # Leukemias and subtypes
  {
    # note, use sample type and anatomic loc vars only, avoid non-cancers from cancer patients
    # leukemia
    ssv <- c("leukemia", "chronic leuk", "chronic myelo",
             "acute leuk", "acute lympho", "acute myel",
             "cml", "CML", "aml", "AML", "ALL")
    mdpost$disease = appendvar("disease", "leukemia", get_filt(get_pstr(ssv), nfilt = T, varl = whichvar))
    mdpost$disease = appendvar("disease", "cancer", get_filt(get_pstr(ssv), nfilt = T, varl = whichvar))
    # chronic leukemias
    ssv <- c("chronic leuk", "chronic leuk", "chronic myelo", "chronic lympho", "cml", "cll", "CML", "CLL")
    mdpost$disease = appendvar("disease", "chronic_leukemia", get_filt(get_pstr(ssv), nfilt = T, varl = whichvar))
    # CML
    ssv <- c("chronic myelo", "cml", "CML")
    mdpost$disease = appendvar("disease", "chronic_myeloid_leukemia", get_filt(get_pstr(ssv), nfilt = T, varl = whichvar))
    # CLL
    ssv <- c("chronic lympho", "cll", "CLL")
    mdpost$disease = appendvar("disease", "chronic_lymphoblastic_leukemia", get_filt(get_pstr(ssv), nfilt = T, varl = whichvar))
    # acute leukemias
    ssv <- c("acute leuk","acute lympho", "acute myel", "aml", "AML", "ALL")
    mdpost$disease = appendvar("disease", "acute_leukemia", get_filt(get_pstr(ssv), nfilt = T, varl = whichvar))
    # AML
    ssv <- c("acute myel", "aml", "AML")
    mdpost$disease = appendvar("disease", "acute_myeloid_leukemia", get_filt(get_pstr(ssv), nfilt = T, varl = whichvar))
    # ALL
    ssv <- c("acute lympho", "ALL")
    mdpost$disease = appendvar("disease", "acute_lymphoblastic_leukemia", get_filt(get_pstr(ssv), nfilt = T, varl = whichvar))
  }
}

save(mdpost, file = nfn)

#-------------------
# tissue annotations
#-------------------
{
  # Tissue position, relation
  {
    # adjacent
    ssv = c("adjacent")
    mdpost$tissue = appendvar("tissue", "adjacent", get_filt(get_pstr(ssv), nfilt = T, varl = whichvar))
    # matched
    ssv = c("match")
    mdpost$tissue = appendvar("tissue", "matched", get_filt(get_pstr(ssv), nfilt = T, varl = whichvar))
    # proximal
    ssv = c("proximal")
    mdpost$tissue = appendvar("tissue", "proximal", get_filt(get_pstr(ssv), nfilt = T, varl = whichvar))
    # distal
    ssv = c("distal")
    mdpost$tissue = appendvar("tissue", "distal", get_filt(get_pstr(ssv), nfilt = T, varl = whichvar))
  }
  
  # GI tract, esophagus to rectum
  {
    # esophageal
    ssv = c("esophag", "oesophag")
    mdpost$tissue = appendvar("tissue", "esophagus", get_filt(get_pstr(ssv)))
    # stomach
    ssv = c("stomach", "gastric")
    mdpost$tissue = appendvar("tissue", "stomach", get_filt(get_pstr(ssv)))
    # small intestine
    ssv = c("small intestine", "small bowel"); filti = get_filt(get_pstr(ssv))
    mdpost$tissue = appendvar("tissue", "small_intestine", filti)
    mdpost$tissue = appendvar("tissue", "intestine", filti)
    # colorectal
    ssv <- c("colorect"); filti = get_filt(get_pstr(ssv))
    mdpost$tissue = appendvar("tissue", "colorectal", filti)
    mdpost$tissue = appendvar("tissue", "intestine", filti)
    # colon
    ssv <- c("colon", "colorec", "large intestine", "cecum"); filti = get_filt(get_pstr(ssv))
    mdpost$tissue = appendvar("tissue", "colon", filti)
    mdpost$tissue = appendvar("tissue", "intestine", filti)
    # rectum
    ssv = c("colorec", "rectal", "rectum", "anus")
    mdpost$tissue = appendvar("tissue", "rectum", get_filt(get_pstr(ssv)))
  }

  # Respiratory system
  {
    # respiratory system
    ssv = c("lung", "bronchi", "alveol", "interstiti", "pleura",
            "trachea", "windpipe", "wind pipe", "bronchi", "airway")
    mdpost$tissue = appendvar("tissue", "respiratory_system", get_filt(get_pstr(ssv)))
    # lung
    ssv = c("lung", "bronchi", "alveol", "interstiti", "pleura")
    mdpost$tissue = appendvar("tissue", "lung", get_filt(get_pstr(ssv)))
    # windpipe
    ssv = c("trachea", "windpipe", "wind pipe", "bronchi", "airway")
    mdpost$tissue = appendvar("tissue", "windpipe", get_filt(get_pstr(ssv)))
  }

  # Nervous system (excluding glial and neuronal categories)
  {
    # respiratory system
    ssv = c("astrocyte", "oligodendrocyte", "ependymal",
            "schwann", "satellite cell")
    mdpost$tissue = appendvar("tissue", "nervous_system", get_filt(get_pstr(ssv)))
  }

  # Chest
  ssv = c("thorax", "chest")
  mdpost$tissue = appendvar("tissue", "chest", get_filt(get_pstr(ssv)))

  # Kidney
  ssv <- c("kidney","renal", "chromophobe", "glomerul", "podocyt", "henle", "glomeruli",
           "nephron", "nephrit", "calyx")
  mdpost$tissue = appendvar("tissue", "kidney", get_filt(get_pstr(ssv)))

  # Liver
  ssv <- c("liver", "hepato", "kupff")
  mdpost$tissue = appendvar("tissue", "liver", get_filt(get_pstr(ssv)))

  # Bladder
  ssv <- c("bladder", "urothel")
  mdpost$tissue = appendvar("tissue", "bladder", get_filt(get_pstr(ssv)))

  # Brain and brain region
  {
    ssv <- c("brain", "cerebel", "dorsolat", "medull", "lobe", "prefront", "occipital",
             "falx", "meningeal", "supratentorial", "fossa", "sellar",
             "grey matter", "gray matter", "white matter")
    mdpost$tissue = appendvar("tissue", "brain", get_filt(get_pstr(ssv)))
    # cerebellum
    ssv <- c("cerebel")
    mdpost$tissue = appendvar("tissue", "cerebellum", get_filt(get_pstr(ssv)))
    # prefrontal lobe
    ssv <- c("prefront")
    mdpost$tissue = appendvar("tissue", "prefrontal_lobe", get_filt(get_pstr(ssv)))
    # grey matter
    ssv <- c("grey matter", "gray matter")
    mdpost$tissue = appendvar("tissue", "gray_matter", get_filt(get_pstr(ssv)))
    # white matter
    ssv <- c("white matter")
    mdpost$tissue = appendvar("tissue", "white_matter", get_filt(get_pstr(ssv)))

  }

  # Embryonic, prenatal tissues
  {
    # placenta
    ssv = c("chorion", "villus", "placent")
    mdpost$tissue = appendvar("tissue", "placenta", get_filt(get_pstr(ssv)))
    # umbilical
    ssv = c("umbilical", "cord")
    mdpost$tissue = appendvar("tissue", "umbilical_cord", get_filt(get_pstr(ssv)))
  }

  # Gametes, gonads, sex-specific, reproductive, inc. cancers
  {
    # uterus
    ssv <- c("uterus", "uteri", "endometr")
    mdpost$tissue = appendvar("tissue", "uterus", get_filt(get_pstr(ssv)))
    # cervix
    ssv <- c("cervix")
    mdpost$tissue = appendvar("tissue", "cervix", get_filt(get_pstr(ssv)))
    # ovary
    ssv <- c("ovary", "ovari")
    mdpost$tissue = appendvar("tissue", "ovary", get_filt(get_pstr(ssv)))
    # vagina
    ssv <- c("vagin")
    mdpost$tissue = appendvar("tissue", "vagina", get_filt(get_pstr(ssv)))
    # labia
    ssv <- c("labia")
    mdpost$tissue = appendvar("tissue", "labia", get_filt(get_pstr(ssv)))
    # fallopian tube
    ssv <- c("fallop")
    mdpost$tissue = appendvar("tissue", "fallopian_tube", get_filt(get_pstr(ssv)))

    # penis
    ssv <- c("penis", "penile")
    mdpost$tissue = appendvar("tissue", "penis", get_filt(get_pstr(ssv)))
    # scrotum
    ssv <- c("scrotum")
    mdpost$tissue = appendvar("tissue", "scrotum", get_filt(get_pstr(ssv)))
    # epididymus
    ssv <- c("epididym")
    mdpost$tissue = appendvar("tissue", "epididymis", get_filt(get_pstr(ssv)))
    # vas deferens
    ssv <- c("vas deferens")
    mdpost$tissue = appendvar("tissue", "vas_deferens", get_filt(get_pstr(ssv)))
    # seminal vesicle
    ssv <- c("seminal vesicle")
    mdpost$tissue = appendvar("tissue", "seminal_vesicle", get_filt(get_pstr(ssv)))
    # prostate
    ssv <- c("prostate")
    mdpost$tissue = appendvar("tissue", "prostate", get_filt(get_pstr(ssv)))
    # testes
    ssv <- c("testic", "teste")
    mdpost$tissue = appendvar("tissue", "testic", get_filt(get_pstr(ssv)))

    # urethra
    ssv <- c("urethra")
    mdpost$tissue = appendvar("tissue", "urethra", get_filt(get_pstr(ssv)))

    # egg
    ssv <- c("egg")
    mdpost$tissue = appendvar("tissue", "egg", get_filt(get_pstr(ssv)))
    # sperm
    ssv <- c("sperm", "spermat")
    mdpost$tissue = appendvar("tissue", "sperm", get_filt(get_pstr(ssv)))
  }

  # Neck, inc. thyroid
  {
    # neck
    ssv = c("neck", "thyroid")
    mdpost$tissue = appendvar("tissue", "neck", get_filt(get_pstr(ssv)))
    # thyroid gland
    ssv = c("thyroid")
    mdpost$tissue = appendvar("tissue", "thyroid_gland", get_filt(get_pstr(ssv)))
  }

  # Eye/optical
  {
    # eye
    ssv = c("eye", "uvea", "optic nerve", "cone", "rod", "retin")
    mdpost$tissue = appendvar("tissue", "eye", get_filt(get_pstr(ssv)))
    # optic nerve
    ssv = c("optic nerve")
    mdpost$tissue = appendvar("tissue", "optic_nerve", get_filt(get_pstr(ssv)))
    # cone
    ssv = c("cone")
    mdpost$tissue = appendvar("tissue", "cone", get_filt(get_pstr(ssv)))
    # rod
    ssv = c("rod")
    mdpost$tissue = appendvar("tissue", "rod", get_filt(get_pstr(ssv)))
    # cancers
    ssv = c("uveal cancer", "uveal melanoma", "retinoblastoma")
    mdpost$tissue = appendvar("tissue", "cancer", get_filt(get_pstr(ssv)))
  }

  # endocrine glands, inc. pancreas and cancers
  {
    # endocrine
    ssv = c("endocrine", "pineal", "pituitary", "pancreas", "pancreat", "adren", "thyroid", "hypothalamus",
            "adrenal cortex", "adreno", "paraganglioma", "paraganglioma", "pheochromocytoma",
            "zona", "glomerulosa", "fasciculata", "reticularis",
            "ovary", "ovari", "testic", "teste")
    mdpost$tissue = appendvar("tissue", "endocrine_system", get_filt(get_pstr(ssv)))
    # pancreas
    ssv = c("pancreas", "pancreat")
    mdpost$tissue = appendvar("tissue", "pancreas", get_filt(get_pstr(ssv)))
  }

  # Skin
  {
    # skin
    ssv <- c("skin", "epidermis", "keratinocyt")
    mdpost$tissue = appendvar("tissue", "skin", get_filt(get_pstr(ssv)))
    # keratinocyte
    ssv <- c("keratinocyt")
    mdpost$tissue = appendvar("tissue", "keratinocyte", get_filt(get_pstr(ssv)))
  }

  # Breast and adipose
  {
    ssv = c("breast")
    mdpost$tissue = appendvar("tissue", "breast", get_filt(get_pstr(ssv)))

    # adipose
    ssv = c("adip", "fat")
    mdpost$tissue = appendvar("tissue", "adipose", get_filt(get_pstr(ssv)))
  }

  # Lymphatic, inc. spleen and thymus
  {
    # lymphatic
    ssv = c("lymph", "spleen", "thymus")
    mdpost$tissue = appendvar("tissue", "lymphatic_system", get_filt(get_pstr(ssv)))
    # thymus
    ssv = c("thymus")
    mdpost$tissue = appendvar("tissue", "thymus", get_filt(get_pstr(ssv)))
    # spleen
    ssv = c("spleen")
    mdpost$tissue = appendvar("tissue", "spleen", get_filt(get_pstr(ssv)))
  }

  # Blood, inc. primary cells
  {
    # blood
    ssv = c("blood", "hematopoiet", "haematopoiet",
            "lymphoid", "myeloid", "natural killer", "nk", "NK",
            "erythrocyte", "mast", "myeloblast", "plasma",
            "monocyte", "lymphocyte", "eosinophil", "neutrophil", "basophil",
            "macrophage", "megakaryocyte", "thrombocyte",
            "wbc", "WBC", "rbc", "RBC",
            "bcell", "b cell", "tcell", "t cell",
            "cd4", "cd5", "cd8", "cytotoxic", "helper")
    mdpost$tissue = appendvar("tissue", "blood", get_filt(get_pstr(ssv)))
    # whole blood
    ssv = c("whole blood")
    mdpost$tissue = appendvar("tissue", "whole_blood", get_filt(get_pstr(ssv)))

    # peripheral blood
    ssv = c("peripheral blood")
    mdpost$tissue = appendvar("tissue", "peripheral_blood", get_filt(get_pstr(ssv)))

    # cord blood
    ssv = c("cord blood")
    mdpost$tissue = appendvar("tissue", "cord_blood", get_filt(get_pstr(ssv)))

    # blood spot
    ssv = c("blood spot")
    mdpost$tissue = appendvar("tissue", "blood_spot", get_filt(get_pstr(ssv)))

    # white blood cells
    ssv = c("wbc", "WBC", "white blood cell",
            "monocyte", "lymphocyte", "eosinophil", "neutrophil", "basophil",
            "bcell", "b cell", "tcell", "t cell",
            "cd4", "cd5", "cd8", "cytotoxic", "helper")
    mdpost$tissue = appendvar("tissue", "white_blood_cell", get_filt(get_pstr(ssv)))

    # Tcells
    ssv = c("tcell", "t cell", "cd4", "cd5", "cd8", "cytotoxic", "helper")
    mdpost$tissue = appendvar("tissue", "t_cell", get_filt(get_pstr(ssv)))
  }

  # Oral, inc. throat
  {
    # oral
    ssv = c("oral", "buccal", "labial", "tongue", "lingual", "throat", "masticatory")
    mdpost$tissue = appendvar("tissue", "oral", get_filt(get_pstr(ssv)))
    # buccal
    ssv = c("buccal")
    mdpost$tissue = appendvar("tissue", "buccal", get_filt(get_pstr(ssv)))
    # throat
    ssv = c("buccal")
    mdpost$tissue = appendvar("tissue", "throat", get_filt(get_pstr(ssv)))
    # tongue
    ssv = c("tongue")
    mdpost$tissue = appendvar("tissue", "tongue", get_filt(get_pstr(ssv)))
  }

  # Nasal
  ssv = c("nasal", "septum")
  mdpost$tissue = appendvar("tissue", "nasal", get_filt(get_pstr(ssv)))

  # Cell types, tissue layers, organ substructures, etc.
  {
    # neurons
    ssv <- c("neur", "nerve", "dendrite", "axon")
    mdpost$tissue = appendvar("tissue", "neuron", get_filt(get_pstr(ssv)))
    # glial cells
    ssv <- c("glia")
    mdpost$tissue = appendvar("tissue", "glia", get_filt(get_pstr(ssv)))
    # epithelial
    ssv = c("epithel")
    mdpost$tissue = appendvar("tissue", "epithelial", get_filt(get_pstr(ssv)))
    # endothelium
    ssv = c("endothel")
    mdpost$tissue = appendvar("tissue", "endothelium", get_filt(get_pstr(ssv)))
    # stem cells
    ssv = c("stem cell", "pluripot", "ipsc", "iPSC")
    mdpost$tissue = appendvar("tissue", "stem_cell", get_filt(get_pstr(ssv)))
    # fibroblast
    ssv = c("fibroblast")
    mdpost$tissue = appendvar("tissue", "fibroblast", get_filt(get_pstr(ssv)))
    # primary cells
    ssv = c("primary cells")
    mdpost$tissue = appendvar("tissue", "primary_cells", get_filt(get_pstr(ssv)))
    # crypt
    ssv = c("crypt")
    mdpost$tissue = appendvar("tissue", "crypt", get_filt(get_pstr(ssv)))
    # muscularis
    ssv = c("muscularis")
    mdpost$tissue = appendvar("tissue", "muscularis", get_filt(get_pstr(ssv)))
    # lamina propria
    ssv = c("lamina propria", "lamina_propria")
    mdpost$tissue = appendvar("tissue", "lamina_propria", get_filt(get_pstr(ssv)))
    # squamous
    ssv = c("squamous")
    mdpost$tissue = appendvar("tissue", "squamous", get_filt(get_pstr(ssv)))
    # ectoderm
    ssv = c("ectoderm")
    mdpost$tissue = appendvar("tissue", "ectoderm", get_filt(get_pstr(ssv)))
    # endoderm
    ssv = c("endoderm")
    mdpost$tissue = appendvar("tissue", "endoderm", get_filt(get_pstr(ssv)))
    # mesoderm
    ssv = c("mesoderm")
    mdpost$tissue = appendvar("tissue", "mesoderm", get_filt(get_pstr(ssv)))
    # melanocyte
    ssv = c("melanocyte")
    mdpost$tissue = appendvar("tissue", "melanocyte", get_filt(get_pstr(ssv)))
    # mucosa
    ssv = c("mucosa")
    mdpost$tissue = appendvar("tissue", "mucosa", get_filt(get_pstr(ssv)))
    # subcutaneous
    ssv = c("subcutaneous")
    mdpost$tissue = appendvar("tissue", "subcutaneous", get_filt(get_pstr(ssv)))
    # visceral
    ssv = c("visceral")
    mdpost$tissue = appendvar("tissue", "visceral", get_filt(get_pstr(ssv)))
  }
}

save(mdpost, file = nfn)

#-----------------------------------
# sample type, metasra-pipeline pred
#-----------------------------------
# append high-confidence meta-sra pipeline sample type predictions
{
  gsmid = mdmap$gsmid
  stype = gsub("'", "",
               gsub(";.*", "",
                    gsub("^.*'sample type':", "", mdmap$msrap_flatjson)))
  spred = as.numeric(gsub("'", "",
                          gsub(";.*", "",
                               gsub("^.*'sample-type confidence':", "", mdmap$msrap_flatjson))))
  dfmd = data.frame(gsm = gsmid,
                    type = stype,
                    pred = round(as.numeric(spred), digits = 3),
                    stringsAsFactors = F)
  dfmd = dfmd[dfmd$gsm %in% mdpost$gsm,]
  dfmd = rbind(dfmd, data.frame(gsm = mdpost[!mdpost$gsm %in% dfmd$gsm,]$gsm,
                                type = rep("NA", nrow(mdpost[!mdpost$gsm %in% dfmd$gsm,])),
                                pred = rep("NA", nrow(mdpost[!mdpost$gsm %in% dfmd$gsm,])),
                                stringsAsFactors = F))
  dfmd = dfmd[order(match(dfmd$gsm, mdpost$gsm)),]
  identical(dfmd$gsm, mdpost$gsm)

  mdpost$sampletype = paste(paste("msraptype", dfmd$type, sep  =":"),
                            paste("msrapconf", dfmd$pred, sep  =":"),
                            sep = ";")
}

save(mdpost, file = nfn)

#----------------------------------------
# sample type, cellosaurus cell line info
#----------------------------------------
# Note: references files described/generated from 'md_cell-lines.R'
# Note: include the mined cell line name with tag "namem:"
{
  # get the most common cell line types from 'misc'
  mdst = rep("NA", nrow(md))
  miscdat = md$misc;
  cll = ifelse(grepl(".*cell_line.*", miscdat) & !grepl(".*cell_line:NA.*", miscdat),
               gsub("(;|$).*", "", gsub("(^|;).*cell_line: ", "", miscdat)),
                    "NA")
  cll = cll[!cll=="NA"]
  cllf = cll[cll %in% ccf$ID]
  # length(cllf) # [1] 1062
  # filt cellosaurus
  ccff = ccf[ccf$ID %in% cllf,]
  ccff$CA = tolower(substr(ccff$CA, 4, nchar(ccff$CA))) # fix group
  # assign info
  for(r in 1:nrow(ccff)){
    dati = as.character(ccff[r,])
    cf = get_filt(get_pstr(dati[1]))
    if(length(cf[cf]) > 0){
      mdpost$sampletype = appendvar("sampletype", "cell_line", cf)
      mdpost$sampletype = appendvar("sampletype", paste0("ccid:", dati[1]), cf)
      mdpost$sampletype = appendvar("sampletype", paste0("ccacc:", dati[2]), cf)
      mdpost$sampletype = appendvar("sampletype", paste0("cccat:", dati[5]), cf)
    }
    message(r)
  }
}

save(mdpost, file = nfn)

#----------------
# age, inc. units
#----------------
{
  # mine age info from mdpre
  xt = as.data.frame(table(unlist(strsplit(md$misc,";")))); xt = xt[rev(order(xt[,2])),]
  aiterms = xt[grepl(".*age_info.*", xt[,1]), 1]
  agedat = md$age
  # format mined age
  {
    af = md$age # original value
    aqvar = "age_info"
    af = gsub("\\..*", "", gsub(" ", "", gsub(aqvar, "", af))) # rm units, spaces, decimals
    af = ifelse(nchar(af) > 2 | nchar(af) == 0, "NA", af) # filter invalid entries
  }
  # age units
  {
    miscdat = md$misc; mdst = rep("", nrow(md))
    # get mined ids and filt na
    whichaa = which(grepl(".*age_info.*", miscdat))
    # filt terms
    ayr = c("year", "yr", "y")
    ady = "day"
    amo = c("month")
    awk = c("week", "wk")
    aqvar = get_pstr(c(ayr, ady, amo, awk))
    aqyr = get_pstr(ayr); aqdy = get_pstr(ady); aqmo = get_pstr(amo); aqwk = get_pstr(awk)
    for(i in whichaa){
      aii = "NA"
      aui = "NA"
      sai = gsub(";.*", "", gsub(".*(|^;)age_info:", "", miscdat[i]))
      if(grepl(aqvar, sai)){
        auval = ifelse(grepl(aqyr, sai), "years",
                       ifelse(grepl(aqdy, sai), "days",
                              ifelse(grepl(aqmo, sai), "months",
                                     ifelse(grepl(aqwk, sai), "weeks", "NA"))))

      } else{
        sai = agedat[i]
        if(grepl(aqvar, sai)){
          auval = ifelse(grepl(aqyr, sai), "years",
                         ifelse(grepl(aqdy, sai), "days",
                                ifelse(grepl(aqmo, sai), "months",
                                       ifelse(grepl(aqwk, sai), "weeks", "NA"))))
        }
      }
      aui = paste0("unitm:", auval)
      mdst[i] = aui
      message(i)
    }
  }
  # add mined info as 'infom'
  {
    vstr = 'infom:'
    iad = "adult"; ife = c("fetal", "foetal"); ipe = "pediatric"; ipr = "prepubescent"
    iadq = get_pstr(iad); ifeq = get_pstr(ife); ipeq = get_pstr(ipe); iprq = get_pstr(ipr)
    infoq = get_pstr(c(iad, ife, ipe, ipr))
    whichinfo = which(grepl(infoq, agedat)) # 1068
    mdsi = rep("", length(agedat))
    for(i in whichinfo){
      di = agedat[i]
      mdsi[i] = ifelse(grepl(iadq, di), paste0(vstr, "adult"),
                       ifelse(grepl(ifeq, di), paste0(vstr, "fetal"),
                              ifelse(grepl(ipeq, di), paste0(vstr, "pediatric"),
                                     ifelse(grepl(ipr, di), paste0(vstr, "prepubescent"), ""))))
    }
  }
  # export for predage inference
  {
    mdage = mdpost
    mdage$age = paste0("valm:", gsub(" ", "", af))
    mdage$age = ifelse(!mdst=="", paste(mdage$age, mdst, sep = ";"), mdage$age)
    mdage$age = ifelse(!mdsi == "", paste(mdage$age, mdsi, sep = ";"), mdage$age)
    mdage$predage = md$predage
  }
  save(mdage, file = "mdage.rda")
}

# make final age var
# notes on tags:
#   unitm: mined unit
#   infom: mined info
#   valm: mined value
#   unitp: ref-based predicted unit from predage/horvath est.

mdpost$age = mdage$age
mdpost$predage = round(as.numeric(as.character(md$predage)), digits = 2)

save(mdpost, file = nfn)

#----------------
# Sex and predsex
#----------------
{
  mdpost$sex = "NA"
  mdpost$disease = appendvar("disease",
                             "chromosomal_abnormality",
                             get_filt(get_pstr(c("XXX", "XXY")),
                                      varl = "gender"))
  ssv = c("klinefelter", "XXY")
  sfilt = get_pstr(ssv); dfilt = get_filt(sfilt, varl = "gender", nfilt = T)
  mdpost$disease = appendvar("disease", "klinefelter_syndrome", dfilt)
  mdpost$sex = ifelse(grepl(get_pstr(c("female", "f", "FEMALE")), md$gender), "F",
                      ifelse(grepl(get_pstr(c("male", "MALE", "m")), md$gender),
                             "M", "NA"))
  mdpost$predsex = md$predsex
}

save(mdpost, file = nfn)

#----------------
# Cell comp. pred
#----------------
cncellcomp = colnames(md)[grepl(".*predcell.*", colnames(md))]
for(clp in cncellcomp){
  cn = colnames(mdpost)
  mdpost = cbind(mdpost, round(as.numeric(md[,clp]), digits = 2))
  colnames(mdpost) = c(cn, clp)
}

save(mdpost, file = nfn)

#--------------------------------------
# Preparation, inc. Fresh Frozen, FFPE
#--------------------------------------
# from prepd object
{
  mdpost$storage = "NA"
  prepd$storage = prepd$preparation
  prepd$storage = ifelse(!grepl(".*FFPE.*", prepd$storage),
                         "F;frozen", 
                         "FFPE;formalin_fixed_paraffin_embedded")
  prepd = prepd[,c("gsm", "storage")]
  mf = mdpost[,c("gsm", "storage")]
  prepd = rbind(prepd, mf[!mf$gsm %in% prepd$gsm,])
  prepd = prepd[order(match(prepd$gsm, mf$gsm)),]
  identical(prepd$gsm, mf$gsm)
  mdpost$storage = prepd$storage
}

save(mdpost, file = nfn)

#-----
# Save
#-----
save(mdpost, file = nfn)
write.csv(mdpost, file = paste0(substr(nfn, 1, nchar(nfn)-4), ".csv"))
